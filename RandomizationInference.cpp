#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> OLS_RcppArma(std::vector<double> ys,
                                 NumericMatrix Xr) {

  // Adapted from RcppArmadillo package fastLm() function
  // http://dirk.eddelbuettel.com/code/rcpp.armadillo.html
  Rcpp::NumericVector yr = wrap(ys);
  int n = Xr.nrow(), k = Xr.ncol();

  arma::mat X(Xr.begin(), n, k, false); // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);

  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  arma::colvec resid = y - X*coef;            // residuals

  std::vector<double> results(1);
  // SSR
  double ssr = arma::as_scalar( arma::trans(resid)*resid );
  ssr = std::max(ssr, std::numeric_limits<double>::epsilon()*1e1);
  results[0] = ssr;

  return results;
}

// [[Rcpp::export]]
std::vector<double> CausalModel_FindTestStats(
    std::vector<double> y, std::vector<int> x, 
    std::vector<double> z, std::vector<double> t, 
    std::vector<double> neighs,
    std::vector<double> age,
    std::vector<double> gender, 
    std::vector<double> priorrel, 
    int n, int m,
    double delta, double tau, int model){

  // Note that y=y(0) uniformity trial outcomes!! -------------------
  std::vector<double> x_double(x.begin(),x.end());
  std::vector<double> ts;
  
  // For residual sum of squares from linear regression
  NumericMatrix Xssr(n, 7);
  // // Vectorized 
  Rcpp::NumericVector intercept(n,1.0);
  Rcpp::NumericVector x_double_r = wrap(x_double);
  Rcpp::NumericVector z_r = wrap(z);
  Rcpp::NumericVector neighs_r = wrap(neighs);
  Rcpp::NumericVector age_r = wrap(age);
  Rcpp::NumericVector gender_r = wrap(gender);
  Rcpp::NumericVector priorrel_r = wrap(priorrel);
  
  Xssr(_,0) = intercept;
  // treatments
  Xssr(_,1) = x_double_r;
  Xssr(_,2) = z_r;
  // covariates
  Xssr(_,3) = neighs_r;
  Xssr(_,4) = age_r;
  Xssr(_,5) = gender_r;
  Xssr(_,6) = priorrel_r;
  
  std::vector<double> ols_ts = OLS_RcppArma(y, Xssr);
  ts.push_back( 1.0/ols_ts[0] );
  
  return ts;
}

// [[Rcpp::export]]
std::vector<double> CausalModel_transformer(
    std::vector<double> y, std::vector<double> x_double,
    std::vector<double> z, std::vector<double> t, 
    double delta, double tau, int n, int model) {
  
  std::vector<double> y_out(y.begin(),y.end());
  for (int i=0; i<n; ++i) {
    // Different causal models that determine y(x)
    double diff_i = 0.0;
    if (model == 0) { 
      // social influence effect is additive in proportion treated 
      diff_i += (delta*x_double[i])+(tau*z[i]);
    }
    y_out[i] -= diff_i; // returns uniformity trial potential outcome
  }
  return y_out;
}

// [[Rcpp::export]]
List OneCausalModel_edgelist_MC(
    std::vector<double> obs_Ys,
    std::vector<int> obs_Xs,
    std::vector<double> age,
    std::vector<double> gender, 
    std::vector<double> priorrel, 
    IntegerMatrix mc_Xs, // N x |Omega| treatment assignments
    IntegerMatrix interfere_edges, // E x 2 matrix for A_[ij] adjacencies
    std::vector<double> d_H0, // hypothesized values of delta, tau
    std::vector<double> t_H0,
    int model) {

  int n=obs_Xs.size();
  int m=std::accumulate(obs_Xs.begin(),obs_Xs.end(),0);
  int numzs = mc_Xs.ncol(); // number of X vectors
  std::vector<double> obs_X_double(obs_Xs.begin(),obs_Xs.end());
  
  std::vector<std::vector<int> > interfere_sets(n, std::vector<int>(n,0));
  for (int ei=0; ei<interfere_edges.nrow(); ei++) {
    interfere_sets[interfere_edges(ei,0)-1][interfere_edges(ei,1)-1] = 1;
  }

  // Matrix multiplications for number (T) and proportion (Z) assigned to X=1
  std::vector<double> obs_Zs(n, 0.0),obs_Ts(n, 0.0),neigh_size(n,0.0);
  std::vector<std::vector<double> > mc_Zs(numzs, std::vector<double>(n,0.0));
  std::vector<std::vector<double> > mc_Ts(numzs, std::vector<double>(n,0.0));
  std::vector<std::vector<int> > mc_X_(numzs, std::vector<int>(n,0));
  for (int mc=0; mc<numzs; ++mc) {
    for (int i=0; i<n; ++i) {
      mc_X_[mc][i] = mc_Xs(i,mc);
    }
  }
  for (int ii=0; ii<n; ++ii) {
    std::vector<int> Ai = interfere_sets[ii];
    int Ai_sum = std::accumulate(Ai.begin(), Ai.end(), 0);
    neigh_size[ii] = (double) Ai_sum;
    obs_Ts[ii] = (double) std::inner_product(
      Ai.begin(), Ai.end(), obs_Xs.begin(), 0);
    if (Ai_sum>0) {
      obs_Zs[ii] = obs_Ts[ii] / neigh_size[ii];
    }
    for (int imc=0; imc<numzs; ++imc) {
      std::vector<int> mc_Z = mc_X_[imc];
      mc_Ts[imc][ii] = (double) std::inner_product(
        Ai.begin(), Ai.end(), mc_Z.begin(), 0);
      if (Ai_sum>0) {
        mc_Zs[imc][ii] = mc_Ts[imc][ii] / neigh_size[ii];
      }
    }
  }

  // List of lists to store results
  List results;

  for (int d=0; d<d_H0.size(); ++d) {
    for (int t=0; t<t_H0.size(); ++t) {
      // fix each value of delta,tau
      double delta = d_H0[d];
      double tau = t_H0[t];

      // Observed test stats that depend on fixed delta,tau
      std::vector<double> unif_Y0_deltatau = CausalModel_transformer(
        obs_Ys, obs_X_double, obs_Zs, obs_Ts, delta, tau, n, model);
      std::vector<double> obs_ts = CausalModel_FindTestStats(
        unif_Y0_deltatau, obs_Xs, obs_Zs, obs_Ts, neigh_size, 
        age, gender, priorrel, 
        n, m, delta, tau, model);
      int num_ts = obs_ts.size();
      std::vector<double> pv(num_ts, 0.0);

      for (int j=0; j<numzs; ++j) {
        // initialize vectors
        std::vector<double> mc_Z = mc_Zs[j];
        std::vector<double> mc_T = mc_Ts[j];
        std::vector<int> mc_X = mc_X_[j];
        std::vector<double> mc_X_double(mc_X.begin(),mc_X.end());

        // Calculate test statistics for y(0) and z -----------------
        std::vector<double> ts_mc_j = CausalModel_FindTestStats(
          unif_Y0_deltatau, mc_X, mc_Z, mc_T, neigh_size, 
          age, gender, priorrel, 
          n, m, delta, tau, model);

        for (int ts=0; ts<num_ts; ++ts) {
          // contribution to respective p-value
          if (std::abs(ts_mc_j[ts]) >= std::abs(obs_ts[ts])*
              exp(-std::numeric_limits<double>::epsilon()*1e2)) {
            pv[ts] += 1.0/((double) numzs);
          }
        }
      }
      results.push_back( List::create(model,delta,tau,pv) );
    }
  }
  return results;
}