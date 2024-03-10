// Author: Zheng Li
// Date: 2022-07-11
// VINTAGE

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
      vec n; // sample sizes
      int p; // number of cis-SNPs
      int k; // number of eigenvalues used for estimation
      vec z1; // transformed and PVE-adjusted z-scores in the eQTL study
      vec z2; // transformed and PVE-adjusted z-scores in the GWAS
      mat Q; // eigenvectors of the LD matrix
      vec lambdas1; // eigenvalues of the LD matrix in the eQTL study
      vec lambdas2; // eigenvalues of the LD matrix in the GWAS
    } dat;
    
    struct Paras
    {
      mat D; // cov[(beta_1j, beta_2j)']
      mat init_D; // initial values of D
      vec sigma_sq; // residual variances
    } paras;
    
    struct PX_EM_Alg
    {
      vec mu_beta1;
      vec mu_beta2;
      vec S1, S2, S4; // block elements of Sigma_beta
      vec alpha; // expansion parameters 
    } EM;
  
    struct Control
    {
      int maxIterEM;
      double tolEM;
      int save_profile;
      int max_profile;
      bool fail;
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
      const vec &lambdas1, const vec &lambdas2, const vec &n, 
      const int p, const int k)
    {
      dat.Q = Q.cols(0, k-1);
      dat.lambdas1 = lambdas1.subvec(0, k-1);
      dat.lambdas2 = lambdas2.subvec(0, k-1);
      dat.z1 = dat.Q.t() * z1;
      dat.z2 = dat.Q.t() * z2;
      dat.p = p;
      dat.k = k;
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
      control.fail = false;
      control.scenario = scenario;
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
        double detD_inv = 1.0 / (paras.D(0,0) * paras.D(1,1) - 
          paras.D(0,1) * paras.D(1,0));
        vec A1 = dat.n(0) / paras.sigma_sq(0) * dat.lambdas1 + 
          paras.D(1,1) * detD_inv;
        vec A4 = dat.n(1) / paras.sigma_sq(1) * dat.lambdas2 + 
          paras.D(0,0) * detD_inv;
        double A2 = -paras.D(0,1) * detD_inv;
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
          1.0 / paras.D(0,0));
        EM.mu_beta1 = EM.S1 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0);
      }
      else if(control.scenario == "r0")
      {
        EM.S1 = 1.0 / (dat.n(0) / paras.sigma_sq(0) * dat.lambdas1 + 
          1.0 / paras.D(0,0));
        EM.mu_beta1 = EM.S1 % dat.z1 * sqrt(dat.n(0)) / paras.sigma_sq(0);
        EM.S4 = 1.0 / (dat.n(1) / paras.sigma_sq(1) * dat.lambdas2 + 
          1.0 / paras.D(1,1));
        EM.mu_beta2 = EM.S4 % dat.z2 * sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
    }
    
    void M_step()
    {
      // 1.update D
      paras.D(0,0) = accu(EM.mu_beta1 % EM.mu_beta1) + accu(EM.S1);
      if(control.scenario == "r0" || control.scenario == "alt")
      {
        paras.D(1,1) = accu(EM.mu_beta2 % EM.mu_beta2) + accu(EM.S4);
      }
      if(control.scenario == "alt")
      {
        paras.D(0,1) = accu(EM.mu_beta1 % EM.mu_beta2) + accu(EM.S2);
        paras.D(1,0) = paras.D(0,1); 
      }
      paras.D = paras.D / dat.k;

      // 2.update expansion parameters alpha1 and alpha2
      EM.alpha(0) = accu(EM.mu_beta1 % dat.z1) / sqrt(dat.n(0));
      EM.alpha(0) /= accu(dat.lambdas1 % EM.mu_beta1 % EM.mu_beta1) +
        accu(dat.lambdas1 % EM.S1);
      if(control.scenario == "alt" || control.scenario == "r0")
      {
        EM.alpha(1) = accu(EM.mu_beta2 % dat.z2) / sqrt(dat.n(1));
        EM.alpha(1) /= accu(dat.lambdas2 % EM.mu_beta2 % EM.mu_beta2) +
          accu(dat.lambdas2 % EM.S4);
      }
      
      // 3.update sigma1_sq and sigma2_sq
      paras.sigma_sq(0) = 1.0 + EM.alpha(0) * EM.alpha(0) * accu(EM.mu_beta1 % 
        dat.lambdas1 % EM.mu_beta1) - 2.0 * EM.alpha(0) * accu(EM.mu_beta1 % 
        dat.z1) / sqrt(dat.n(0));
      paras.sigma_sq(0) += EM.alpha(0) * EM.alpha(0) * 
        accu(EM.S1 % dat.lambdas1);
      if(control.scenario == "alt" || control.scenario == "r0")
      {
        paras.sigma_sq(1) = 1.0 + EM.alpha(1) * EM.alpha(1) * accu(EM.mu_beta2 % 
          dat.lambdas2 % EM.mu_beta2) - 2.0 * EM.alpha(1) * accu(EM.mu_beta2 % 
          dat.z2) / sqrt(dat.n(1));
        paras.sigma_sq(1) += EM.alpha(1) * EM.alpha(1) * 
          accu(EM.S4 % dat.lambdas2);
      }
      
      // handle edge cases
      if((paras.D(0,0) * dat.p > 1.0) || (paras.D(1,1) * dat.p) > 1.0 ||
        paras.sigma_sq(0) < 0.0 || paras.sigma_sq(0) > 1.0 ||
        paras.sigma_sq(1) < 0.0 || paras.sigma_sq(1) > 1.0)
      {
        control.fail = true;
        paras.D(0,0) = 0.01 / dat.p;
        paras.D(1,1) = 0.0004 / dat.p;
        paras.sigma_sq(0) = 1.0 - paras.D(0,0) * dat.k;
        paras.sigma_sq(1) = 1.0 - paras.D(1,1) * dat.k;
        EM.alpha.fill(1.0);
      }
    }
    
    void reduction_step()
    {
      paras.D = diagmat(EM.alpha) * paras.D * diagmat(EM.alpha);
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
        profile.loglik -= dat.k * logdet;
        profile.loglik += accu(log(EM.S4 - EM.S2 % EM.S2 / EM.S1)) + 
          accu(log(EM.S1));
        profile.loglik += accu(EM.mu_beta2 % dat.z2) *
          sqrt(dat.n(1)) / paras.sigma_sq(1);
      }
      else if(control.scenario == "null")
      {
        profile.loglik -= dat.k * log(paras.D(0,0));
        profile.loglik += accu(log(EM.S1));
      }
      else if(control.scenario == "r0")
      {
        profile.loglik -= dat.k * (log(paras.D(0,0)) + log(paras.D(1,1)));
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
          if(control.fail || (abs(profile.logliks(profile_iter-1) - 
            profile.logliks(profile_iter-2)) < control.tolEM))
          {
            break;
          }
        }
        M_step();
        reduction_step();
        iter += 1;
      }
    }
    
    List genetic_var_test()
    {
      vec scores(2), test_stats, pvalues(11);
      mat I(2, 2); // information matrix of theta1 under H0

      // evaluate information
      // note: I and scores have been reconciled to avoid integer overflow
      I(0,0) = accu(dat.lambdas2 % dat.lambdas2) / 2.0;
      I(0,0) = I(0,0) - dat.k * dat.k / dat.n(1) / 2.0;
      I(1,1) = accu(dat.lambdas1 % EM.S1 % dat.lambdas2) / 
        paras.sigma_sq(0) / paras.D(0,0);
      
      // evaluate score
      scores(0) = accu(dat.z2 % dat.z2) / 2.0;
      scores(1) = accu(dat.z1 % EM.S1 % dat.z2) / 
        paras.sigma_sq(0) / paras.D(0,0);
      scores(0) /= sqrt(I(0,0));
      scores(1) /= sqrt(I(1,1));

      // linearly weighted score test statistics
      test_stats = regspace(1, -0.1, 0) * scores(0) + 
        regspace(0, 0.1, 1) * abs(scores(1));
      
      // simulate test statistics
      int siter = 0;
      mat sim_test_stats(11, testing.B);
      vec temp = dat.z1 % EM.S1 / paras.sigma_sq(0) / paras.D(0,0);
      for(int i = 0; i < testing.B; i++)
      {
        vec sim_scores(2);
        vec sim_z2(dat.k, fill::randn);
        sim_z2 = sim_z2 % sqrt(dat.lambdas2);
        
        sim_scores(0) = accu(sim_z2 % sim_z2) / 2.0;
        sim_scores(1) =  accu(temp % sim_z2);
        sim_scores(0) /= sqrt(I(0,0));
        sim_scores(1) /= sqrt(I(1,1));
        sim_test_stats.col(i) = regspace(1, -0.1, 0) * sim_scores(0) 
          + regspace(0, 0.1, 1) * abs(sim_scores(1));
      }
      
      // evaluate p-values using simulated test statistics
      for(int i = 0; i < 11; i++)
      {
        pvalues(i) = (double)sum(test_stats(i) < sim_test_stats.row(i)) / 
          testing.B;
      }
      
      List output = List::create(
        _["I"] = I,
        _["B"] = testing.B,
        _["test_stats"] = test_stats,
        _["pvalues"] = pvalues
      );
      return output;
    }
    
    List genetic_corr_test()
    {
      double u = accu(dat.z1 % EM.S1 % dat.z2 % EM.S4) / 
        prod(paras.sigma_sq) / paras.D(0,0) / paras.D(1,1);
      double var_u = accu(EM.S1 % EM.S4 % dat.lambdas1 % dat.lambdas2) / 
        prod(paras.sigma_sq) / paras.D(0,0) / paras.D(1,1);
      double Tr = u * u / var_u;
      
      List output = List::create(
        _["u"] = u,
        _["var_u"] = var_u,
        _["Tr"] = Tr
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
        _["fail"] = control.fail
      );
      if(control.scenario == "alt")
      {
        r = paras.D(0,1) / sqrt(paras.D(0,0)) / sqrt(paras.D(1,1));
        output.push_front(r, "r");
        output.push_front(paras.D(1,1) * dat.p, "h2_sq");
      }
      else if(control.scenario == "r0")
      {
        output.push_front(paras.D(1,1) * dat.p, "h2_sq");
        List corr_test_res = genetic_corr_test();
        output.push_back(corr_test_res, "corr_test");
      }
      output.push_front(paras.D(0,0) * dat.p, "h1_sq");
      
      if(control.scenario == "null")
      {
        List var_test_res = genetic_var_test();
        output.push_back(var_test_res, "var_test");
      }
      return output;
    }
};


// [[Rcpp::export]]
List vintage(const arma::vec &z1, const arma::vec &z2, const arma::mat Q, 
  const arma::vec lambdas1, const arma::vec lambdas2, const arma::vec &n, 
  const int p, const int k, const arma::mat init_D, const int maxIterEM, 
  const double tolEM, const std::string scenario, const int B, 
  const int save_profile = 1)
{
  wall_clock timer;
  timer.tic();
  VINTAGEModel model;
  
  model.load_data(z1, z2, Q, lambdas1, lambdas2, n, p, k);
  model.set_control(maxIterEM, tolEM, save_profile, scenario);
  model.set_testing(B);
  model.init_paras(init_D);
  model.run();
  
  List output = model.get_output();
  double elapsed = timer.toc();
  output.push_back(elapsed, "elapsed_time");
  return output;
}

