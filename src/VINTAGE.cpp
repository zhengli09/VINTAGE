// Author: Zheng Li
// Date: 2022-07-11
// Bivariate linear mixed model that intends
// to integrate TWAS and SKAT

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

class VINTAGEModel
{
  private:
    struct Data
    {
      vec n; // n1 and n2
      int p;
      vec z1; // transformed z-scores of eQTL study
      vec z2; // transformed z-scores of GWAS
      mat Q; // eigenvectors of the SNP correlation matrix
      vec lambdas1; // eigenvalues of the SNP correlation matrix in eQTL study
      vec lambdas2; // eigenvalues of the SNP correlation matrix in GWAS
    } dat;
    
    struct Paras
    {
      mat D; // cov[(beta_1j, beta_2j)']=[h1_sq&rho\\rho&h2_sq]/p
      mat init_D; // store initial values of D from LDSC
      vec sigma_sq; // residual variance 
    } paras;
    
    struct PX_EM_Alg
    {
      vec mu_beta1;
      vec mu_beta2;
      vec S1, S2, S4; // block elements of Sigma_beta
      vec alpha; // expanded parameter for faster convergence
    } EM;
  
    struct Control
    {
      int maxIterEM;
      double tolEM;
      int save_profile;
      int max_profile;
      bool is_converged;
      bool is_D_fixed; // fix D to be LDSC estimates if it is not within 0-1
      string scenario;
    } control;
    
    struct Testing
    {
      int B; // number of simulations used for constructing the null
    } testing;
    
    struct Profile
    {
      double loglik;
      vec logliks;
    } profile;
  
    int iter = 0;
    int profile_iter = 0;
  
  public:
    void load_data(const vec &z1, const vec &z2, const mat &Q, 
      const vec &lambdas1, const vec &lambdas2, const vec &n)
    {
      dat.Q = Q;
      dat.lambdas1 = lambdas1;
      dat.lambdas2 = lambdas2;
      dat.z1 = dat.Q.t() * z1;
      dat.z2 = dat.Q.t() * z2;
      dat.p = dat.z1.n_elem;
      dat.n = n;
    }
    
    void set_control(const int maxIterEM, const double tolEM, 
      const int save_profile, const string scenario)
    {
      control.maxIterEM = maxIterEM;
      control.tolEM = tolEM;
      control.save_profile = save_profile;
      if(control.save_profile > 0)
      {
        control.max_profile = control.maxIterEM / control.save_profile + 1;
      }
      control.scenario = scenario;
      control.is_converged = false;
      control.is_D_fixed = false;
    }
    
    void set_testing(const int B)
    {
      testing.B = B;
    }
    
    void init_paras(mat D)
    {
      paras.init_D = D;
      paras.D = paras.init_D;
      paras.sigma_sq = vec(2, fill::ones);
      EM.alpha = vec(2, fill::ones);
      profile.logliks.zeros(control.max_profile);
    }
    
    // PX-EM algorithm
    void E_step()
    {
      if(control.scenario == "alt")
      {
        // update Sigma_beta
        double detD = dat.p / (paras.D(0,0) * paras.D(1,1) - 
          paras.D(0,1) * paras.D(1,0));
        vec A1 = dat.n(0) / paras.sigma_sq(0) * dat.lambdas1 + 
          paras.D(1,1) * detD;
        vec A4 = dat.n(1) / paras.sigma_sq(1) * dat.lambdas2 + 
          paras.D(0,0) * detD;
        double A2 = -paras.D(0,1) * detD;
        EM.S1 = 1.0 / (A1 - A2 * A2 / A4);
        EM.S4 = 1.0 / (A4 - A2 * A2 / A1);
        EM.S2 = -A2 * EM.S1 / A4;
        
        // update mu_beta
        EM.mu_beta1 = EM.S1 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0) +
          EM.S2 % dat.z2 * sqrt(dat.n(1)) / paras.sigma_sq(1);
        EM.mu_beta2 = EM.S2 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0) +
          EM.S4 % dat.z2 * sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
      else if(control.scenario == "null")
      {
        EM.S1 = 1.0 / (dat.n(0) / paras.sigma_sq(0) * dat.lambdas1 + 
          dat.p / paras.D(0,0));
        EM.mu_beta1 = EM.S1 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0);
      }
      else if(control.scenario == "r0")
      {
        EM.S1 = 1.0 / (dat.n(0) / paras.sigma_sq(0) * dat.lambdas1 + 
          dat.p / paras.D(0,0));
        EM.mu_beta1 = EM.S1 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0);
        EM.S4 = 1.0 / (dat.n(1) / paras.sigma_sq(1) * dat.lambdas2 + 
          dat.p / paras.D(1,1));
        EM.mu_beta2 = EM.S4 % dat.z2 * sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
    }
    
    void M_step()
    {
      // update D
      if(!control.is_D_fixed)
      {
        paras.D(0,0) = accu(EM.mu_beta1 % EM.mu_beta1) + accu(EM.S1);
        if(control.scenario == "alt")
        {
          paras.D(1,1) = accu(EM.mu_beta2 % EM.mu_beta2) + accu(EM.S4);
          paras.D(0,1) = accu(EM.mu_beta1 % EM.mu_beta2) + accu(EM.S2);
          paras.D(1,0) = paras.D(0,1); 
        }
        else if(control.scenario == "r0")
        {
          paras.D(1,1) = accu(EM.mu_beta2 % EM.mu_beta2) + accu(EM.S4);
        }
      }
      if(paras.D(0,0) > 1.0 || paras.D(1,1) > 1.0)
      {
        cout << "Warning: Algorithm is unstable likely due to the LD mismatch.";
        cout << " Fix the D matrix using LDSC estimates" << endl;
        control.is_D_fixed = true;
        paras.D = paras.init_D;
      }
      
      // update expansion parameters alpha1 and alpha2
      EM.alpha(0) = accu(EM.mu_beta1 % dat.z1) / sqrt(dat.n(0));
      EM.alpha(0) /= accu(dat.lambdas1 % EM.mu_beta1 % EM.mu_beta1) +
        accu(dat.lambdas1 % EM.S1);
      if(control.scenario == "alt" || control.scenario == "r0")
      {
        EM.alpha(1) = accu(EM.mu_beta2 % dat.z2) / sqrt(dat.n(1));
        EM.alpha(1) /= accu(dat.lambdas2 % EM.mu_beta2 % EM.mu_beta2) +
          accu(dat.lambdas2 % EM.S4);
      }
      
      // update sigma1_sq and sigma2_sq
      paras.sigma_sq(0) = 1 + EM.alpha(0) * EM.alpha(0) * accu(EM.mu_beta1 % 
        dat.lambdas1 % EM.mu_beta1) - 2 * EM.alpha(0) * accu(EM.mu_beta1 % 
        dat.z1) / sqrt(dat.n(0));
      paras.sigma_sq(0) = max(0.0, paras.sigma_sq(0));
      paras.sigma_sq(0) += EM.alpha(0) * EM.alpha(0) * 
        accu(EM.S1 % dat.lambdas1);
      if(control.scenario == "alt" || control.scenario == "r0")
      {
        paras.sigma_sq(1) = 1 + EM.alpha(1) * EM.alpha(1) * accu(EM.mu_beta2 % 
          dat.lambdas2 % EM.mu_beta2) - 2.0 * EM.alpha(1) * accu(EM.mu_beta2 % 
          dat.z2) / sqrt(dat.n(1));
        paras.sigma_sq(1) = max(0.0, paras.sigma_sq(1));
        paras.sigma_sq(1) += EM.alpha(1) * EM.alpha(1) * 
          accu(EM.S4 % dat.lambdas2);
      }
      // constrain sigma2_sq to be larger than 0.1
      paras.sigma_sq.transform([](double x){return (x < 0.1 ? 0.1 : x);});
    }
    
    void reduction_step()
    {
      if(!control.is_D_fixed)
      {
        paras.D = diagmat(EM.alpha) * paras.D * diagmat(EM.alpha);
      }
    }
    
    void update_loglik()
    {
      double logdet;
      double sign;
      profile.loglik = -dat.n(0) * log(paras.sigma_sq(0));
      profile.loglik -= dat.n(1) * log(paras.sigma_sq(1));
      profile.loglik -= dat.n(0) / paras.sigma_sq(0);
      profile.loglik -= dat.n(1) / paras.sigma_sq(1);
      profile.loglik += accu(EM.mu_beta1 % dat.z1) * 
        sqrt(dat.n(0)) / paras.sigma_sq(0);
      if(control.scenario == "alt")
      {
        log_det(logdet, sign, paras.D);
        profile.loglik -= dat.p * logdet;
        profile.loglik += accu(log(EM.S4 - EM.S2 % EM.S2 / EM.S1)) + 
          accu(log(EM.S1));
        profile.loglik += accu(EM.mu_beta2 % dat.z2) *
          sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
      else if(control.scenario == "null")
      {
        profile.loglik -= dat.p * log(paras.D(0,0));
        profile.loglik += accu(log(EM.S1));
      }
      else if(control.scenario == "r0")
      {
        profile.loglik -= dat.p * (log(paras.D(0,0)) + log(paras.D(1,1)));
        profile.loglik += accu(log(EM.S1)) + accu(log(EM.S4));
        profile.loglik += accu(EM.mu_beta2 % dat.z2) *
          sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
    }
    
    void save_profile()
    {
      if((iter+1) % control.save_profile == 0)
      {
        update_loglik();
        profile.logliks(profile_iter) = profile.loglik;
        profile_iter += 1;
      }
    }
    
    void monitor()
    {
      if((iter+1) % control.save_profile == 0)
      {
        cout << "Iter " << iter+1 << ": loglik = " << profile.loglik << endl;
      }
    }
    
    void run()
    {
      while(iter < control.maxIterEM)
      {
        E_step();
        save_profile();
        monitor();
        if(((iter+1) % control.save_profile == 0) && profile_iter > 1)
        {
          if(abs(profile.logliks(profile_iter-1) - 
            profile.logliks(profile_iter-2)) < control.tolEM)
          {
            control.is_converged = true;
            break;
          }
        }
        M_step();
        reduction_step();
        iter += 1;
      }
    }
    
    List sim_based_test()
    {
      cout << "Start simulation-based testing..." << endl;
      vec scores(2), test_stats, pvalues(11);
      mat I(2, 2); // information matrix of theta1 under H0
      // components
      vec LLS = dat.lambdas1 % dat.lambdas1 % EM.S1;
      vec LS = dat.lambdas1 % EM.S1;
 
      // evaluate information
      I(0,0) = dat.n(1) * dat.n(1) * accu(dat.lambdas2 % dat.lambdas2) / 2.0;
      I(0,0) = I(0,0) - dat.n(1) * dat.p * dat.p / 2.0;
      I(1,1) = accu(dat.lambdas1 % dat.lambdas2) - accu(LLS % dat.lambdas2) * 
        dat.n(0) / paras.sigma_sq(0);
      I(1,1) = I(1,1) * prod(dat.n) / paras.sigma_sq(0);
      
      // evaluate score
      scores(0) = accu(dat.z2 % dat.z2) * dat.n(1) / 2.0;
      scores(1) = (accu(dat.z1 % dat.z2) - sqrt(dat.n(0)) * accu(EM.mu_beta1 % 
        dat.lambdas1 % dat.z2)) * sqrt(prod(dat.n)) / paras.sigma_sq(0);
      scores(0) /= sqrt(I(0,0));
      scores(1) /= sqrt(I(1,1));

      // linearly weighted score test statistics
      test_stats = regspace(1, -0.1, 0) * scores(0) + 
        regspace(0, 0.1, 1) * abs(scores(1));
      
      // simulate test statistics
      int siter = 0;
      mat sim_test_stats(11, testing.B);
      vec temp = dat.z1 % (1.0 - LS * dat.n(0) / paras.sigma_sq(0)) * 
        sqrt(prod(dat.n)) / paras.sigma_sq(0);
      for(int i = 0; i < testing.B; i++)
      {
        vec sim_scores(2);
        vec sim_z2(dat.p, fill::randn);
        sim_z2 = sim_z2 % sqrt(dat.lambdas2);
        
        sim_scores(0) = accu(sim_z2 % sim_z2) * dat.n(1) / 2.0;
        sim_scores(1) =  accu(temp % sim_z2);
        sim_scores(0) /= sqrt(I(0,0));
        sim_scores(1) /= sqrt(I(1,1));
        sim_test_stats.col(i) = regspace(1, -0.1, 0) * sim_scores(0) 
          + regspace(0, 0.1, 1) * abs(sim_scores(1));
      }
      
      // evaluate p-values using simulated test statistics
      for(int i = 0; i < 11; i++)
      {
        pvalues(i) = ((double)sum(test_stats(i) < sim_test_stats.row(i)) +
          1.0) / (testing.B + 1.0);
      }
      
      List output = List::create(
        _["I"] = I,
        _["B"] = testing.B,
        _["test_stats"] = test_stats,
        _["pvalues"] = pvalues
      );
      return output;
    }
    
    // List test_h2_sq()
    // {
    //   double u1, u2, var_u1, var_u2, w1, w2;
    //   double score_pos, score_neg, e, nu, kappa;
    //   
    //   u1 = accu(dat.z2 % dat.z2) * dat.n(1) / 2.0;
    //   u2 = accu(dat.z1 % dat.z2) - sqrt(dat.n(0)) * accu(EM.mu_beta1 % 
    //     dat.lambdas1 % dat.z2);
    //   u2 = u2 * sqrt(prod(dat.n)) / paras.sigma_sq(0);
    // 
    //   var_u1 = dat.n(1) * dat.n(1) * accu(dat.lambdas2 % dat.lambdas2) / 2.0;
    //   var_u2 = accu(dat.lambdas1 % dat.lambdas2) - accu(dat.lambdas1 % 
    //     dat.lambdas2 % dat.lambdas1 % EM.S1) * dat.n(0) / paras.sigma_sq(0);
    //   var_u2 = var_u2 * prod(dat.n) / paras.sigma_sq(0);
    // 
    //   w1 = 1.0 / sqrt(var_u1);
    //   w2 = 1.0 / sqrt(var_u2);
    //   e = w1 * dat.n(1) * dat.p / 2.0;
    //   nu = e * e;
    //   kappa = 1.0 / e;
    //   score_pos = (w1 * u1 + w2 * u2) / kappa;
    //   score_neg = (w1 * u1 - w2 * u2) / kappa;
    //   
    //   List output = List::create(
    //     _["score_pos"] = score_pos,
    //     _["score_neg"] = score_neg,
    //     _["nu"] = nu,
    //     _["u1"] = u1,
    //     _["u2"] = u2,
    //     _["w1"] = w1,
    //     _["w2"] = w2
    //   );
    //   return output;
    // }
    
    List test_r0()
    {
      double u;
      vec eigvals;
      double var_u;
      u = accu((dat.z1 - sqrt(dat.n(0)) * dat.lambdas1 % EM.mu_beta1) %
        (dat.z2 - sqrt(dat.n(1)) * dat.lambdas2 % EM.mu_beta2));
      u = u * sqrt(prod(dat.n)) / prod(paras.sigma_sq); 
      eigvals = dat.lambdas1 - dat.lambdas1 % EM.S1 % dat.lambdas1 * dat.n(0) /
        paras.sigma_sq(0);
      eigvals %= dat.lambdas2 - dat.lambdas2 % EM.S4 % dat.lambdas2 * dat.n(1) /
        paras.sigma_sq(1);
      eigvals *= prod(dat.n) / prod(paras.sigma_sq);
      var_u = accu(eigvals);
      eigvals = sqrt(eigvals) / 2.0;
      eigvals = join_cols(-eigvals, eigvals);
      
      List output = List::create(
        _["u"] = u,
        _["eigvals"] = eigvals,
        _["var_u"] = var_u
      );
      return output;
    }
    
    List get_output()
    {
      List output;
      double r;
      output = List::create(
        _["sigma1_sq"] = paras.sigma_sq(0),
        _["sigma2_sq"] = paras.sigma_sq(1),
        _["logliks"] = profile.logliks.head(profile_iter),
        _["is_D_fixed"] = control.is_D_fixed,
        _["is_converged"] = control.is_converged
      );
      if(control.scenario == "alt")
      {
        r = paras.D(0,1) / sqrt(paras.D(0,0)) / sqrt(paras.D(1,1));
        output.push_front(r, "r");
        output.push_front(paras.D(1,1), "h2_sq");
      }
      else if(control.scenario == "r0")
      {
        output.push_front(paras.D(1,1), "h2_sq");
        List testing_r0 = test_r0();
        output.push_back(testing_r0, "test_r0");
      }
      output.push_front(paras.D(0,0), "h1_sq");
      
      if(control.scenario == "null")
      {
        // List testing_h2_sq = test_h2_sq();
        // output.push_back(testing_h2_sq, "test_h2");
        List test_h2_sq = sim_based_test();
        output.push_back(test_h2_sq, "test_h2");
      }
      return output;
    }
};


// [[Rcpp::export]]
List vintage(const arma::vec &z1, const arma::vec &z2, const arma::mat Q, 
  const arma::vec lambdas1, const arma::vec lambdas2, const arma::vec &n, 
  const arma::mat init_D, const int maxIterEM, const double tolEM, 
  const std::string scenario, const int B, const int save_profile = 10)
{
  wall_clock timer;
  timer.tic();
  VINTAGEModel model;
  
  model.load_data(z1, z2, Q, lambdas1, lambdas2, n);
  model.set_control(maxIterEM, tolEM, save_profile, scenario);
  model.set_testing(B);
  model.init_paras(init_D);
  model.run();
  
  List output = model.get_output();
  double elapsed = timer.toc();
  output.push_back(elapsed, "elapsed_time");
  return output;
}

