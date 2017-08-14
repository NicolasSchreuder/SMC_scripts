// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <random>
#include <iostream>
#include <R.h>
#include <math.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions/normal.hpp>
using boost::math::normal;

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "sqmc/R_Functions/include/SamplePack/DigitalNetsBase2.h"
#include "sqmc/R_Functions/include/Generate_RQMC.hpp"

  

using namespace Rcpp;

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
arma::vec StratifiedSampling(const double n){
  // auto-cast in argument (from int to double)
  // Stratified sampling over [0,1], returns n sorted uniforms
  std::mt19937 generator(134);
  arma::vec x = arma::zeros(n);
  for(int i=0; i<n; i++){
    double a = i/n;
    double b = (i+1)/n;
    std::uniform_real_distribution<double> unif(a,b);
    x(i) = unif(generator);
  }
  return x;
};

// [[Rcpp::export]]
arma::vec Multinomial(int n, arma::vec weights){
  arma::vec unif = StratifiedSampling(n);
  double sum_weights = weights(0);
  int k = 0;
  arma::vec sample = arma::zeros(n);
  
  for(int i=0; i<n;i++){
    while (unif[i] > sum_weights){
      k++;
      sum_weights += weights[k];
    }
    sample[i] = k;
  }
  return sample;
}

// [[Rcpp::export]]
arma::vec Sample(arma::mat particles_weights, const int nb_particles, const int time){
  // Return multinomial sample
  return Multinomial(nb_particles, particles_weights.col(time));
}

// [[Rcpp::export]]
double G(const double &y, const double &x, const double &sigma){
  // Density of Y | X=x
  normal N(x,sigma);
  return pdf(N,y);
}

// [[Rcpp::export]]
double F(const double &x, const double &mu, const double &sigma){
  // Kernel associated to the markov chain X_t
  normal N(mu, sigma);
  return pdf(N, x);
}


// [[Rcpp::export]]
arma::vec SMC(arma::vec y, const int M, const int T){
  
  set_seed(42);
  boost::mt19937 rng;
  
  arma::mat particles(M, T, arma::fill::zeros);
  arma::mat weights(M, T, arma::fill::zeros);
  
  // Initial particles
  for(int i=0; i<M; i++){
    particles(i,0) = R::rnorm(0, 1);
    weights(i,0) = G(y(0), particles(i,0), 1.)*F(particles(i,0),0.,1.);
  }

  // Normalization
  double sum_weights = sum(weights.col(0));
  for(int i=0; i<M; i++){
    weights(i,0) /= sum_weights;
  }
  
  arma::vec x_pred(T);

  x_pred(0) = dot(weights.col(0), particles.col(0));
  
  for(int t=1; t<T; t++){
    arma::vec resampling = Multinomial(M, weights.col(t-1));
    
    for(int m=0; m<M; m++){
      // boost::normal_distribution<> nd(particles(resampling(m), t-1), 1.0);
      // boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
      // particles(m, t) = var_nor();
      particles(m, t) = R::rnorm(particles(resampling(m), t-1), 1.);
      weights(m,t) = G(y(t), particles(m,t), 1.);
    }
    
    sum_weights = sum(weights.col(t));
    for(int m=0; m<M; m++){
      weights(m,t) /= sum_weights;
    }
    
    x_pred(t) = dot(weights.col(t), particles.col(t));
  }
  
  
  // Backward smoothing/sampling
  arma::mat backward_particles(M, T, arma::fill::zeros);
  arma::vec resampling = Multinomial(M, weights.col(T-1));
  
  for(int m=0; m<M; m++){
    backward_particles(m, T-1) = particles(resampling(m), T-1);
  }
  
  arma::vec wgts(M);
  
  for(int t=T-2; t>=0; t--){
    double sum = 0;
    for(int m=0; m<M; m++){
      wgts(m) = weights(m,t)*F(backward_particles(m,t+1), particles(m,t), 1.);
      sum += wgts(m);
    }
    
    for(int m=0; m<M; m++){
      wgts(m) /= sum;
    }
    
    resampling = Multinomial(M, wgts);
    
    for(int m=0; m<M; m++){
      backward_particles(m, t) = particles(resampling(m), t);
    }
  }
  arma::vec x_smooth(T);
  for(int t=0; t<T; t++){
    x_smooth[t] = mean(backward_particles.col(t));
  }
  return x_smooth;
}


// Filtre interne 
// [[Rcpp::export]]
double SMC_example_1(arma::vec y, const int M, const int D, arma::mat Sigma, arma::vec x_prev, double sigma_y){
  // Input //
  // y : data
  // M : number of particles
  // D : dimension of data
  // Sigma : Variance-covariance matrix
  // x_prev : previous x
  // sigmay_y : variance of data
  
  set_seed(42);
  boost::mt19937 rng;
  boost::normal_distribution<double> normdist(0., 1.);
  
  // Initialize particles and weights
  arma::mat particles(M, D, arma::fill::zeros);
  arma::mat weights(M, D, arma::fill::zeros);
  
  for(int i=0; i < M; i++){
    particles(i, 0) = 0.5 * x_prev(0) + sqrt(Sigma(0, 0)) * normdist(rng);
    weights(i, 0) = G(y(0), particles(i, 0), sigma_y);
  }
  
  double log_likelihood = log(mean(weights.col(0)));
  
  // Normalization
  double sum_weights = sum(weights.col(0));
  for(int i=0; i<M; i++){
    weights(i,0) = weights(i,0)/sum_weights;
  }
  
  for(int d=1; d<D; d++){
    // Resampling
    arma::vec resampling = Multinomial(M, weights.col(d-1));
    
    // Standard deviation 
    double sd = Sigma(d, d) - pow(Sigma(d-1, d), 2)/Sigma(d-1, d-1);
    
    for(int m = 0; m < M; m++){
      double mu = 0.5*x_prev(d) + Sigma(d-1, d)/Sigma(d-1, d-1)*(particles(resampling(m), d-1) - 0.5*x_prev(d-1));
      particles(m, d) = sqrt(sd) * normdist(rng) + mu;
      weights(m,d) = G(y(d), particles(m, d), sigma_y);
    }
    
    log_likelihood += log(mean(weights.col(d)));
    
    sum_weights = sum(weights.col(d));
    
    for(int m = 0; m < M; m++){
      weights(m,d) /= sum_weights;
    }
  }
  
  return log_likelihood;
}

// [[Rcpp::export]]
arma::vec resampling_QMC(int n, arma::vec weights, arma::vec sobol){
  
  double sum_weights = weights(0);
  int k = 0;
  arma::vec sample = arma::zeros(n);
  
  for(int i=0; i < n;i++){
    while (sobol(i) > sum_weights){
      k++;
      sum_weights += weights(k);
    }
    sample(i) = k;
  }
  return sample;
}

// Filtre interne (coordonnées) avec QMC
// [[Rcpp::export]]
double SMC_example_1_QMC(arma::vec y, const int M, const int D, arma::mat Sigma, arma::vec x_prev, double sigma_y){
  // Input //
  // y : data
  // M : number of particles
  // D : dimension of data
  // Sigma : Variance-covariance matrix
  // x_prev : previous x
  // sobol_sequence : sobol sequence of dimension M*D
  // sigmay_y : variance of data
  
  set_seed(42);

  int seq = 1;
  int N = 10;
  int d = 2;
  
  Scrambled* a = Scrambled_Create(seq, N, d);
  
  arma::mat sobol_sequence;
  
  // Initial particles
  arma::mat particles(M, D, arma::fill::zeros);
  arma::mat weights(M, D, arma::fill::zeros);
  
  boost::math::normal dist(0.0, 1.0);
  
  for(int m=0; m < M; m++){
    particles(m,0) = 0.5*x_prev(0) + sqrt(Sigma(0,0)) * quantile(dist, sobol_sequence(m, 1));
  }
  
  arma::vec particles_sorted = sort(particles.col(0));
  
  for(int m=0; m<M; m++){
    weights(m,0) = G(y(0), particles_sorted(m), sigma_y);
  }
  
  double log_likelihood = log(mean(weights.col(0)));
  
  // Normalization
  double sum_weights = sum(weights.col(0));
  for(int i=0; i<M; i++){
    weights(i,0) /= sum_weights;
  }
  
  for(int d=1; d < D; d++){
    // Resampling
    arma::vec resampling = resampling_QMC(M, weights.col(d-1), sobol_sequence.col(0));

    // Standard deviation 
    double sd = Sigma(d, d) - pow(Sigma(d-1, d), 2)/Sigma(d-1, d-1);
    
    arma::vec sobol_gaussian  = shuffle(sobol_sequence.col(1));
    
    for(int m = 0; m < M; m++){
      double mu = 0.5*x_prev(d) + Sigma(d-1,d)/Sigma(d-1,d-1)*(particles(resampling(m), d-1) - 0.5*x_prev(d-1));
      particles(m, d) = sqrt(sd) * quantile(dist, sobol_gaussian(m)) + mu;
    }
    
    arma::vec particles_sorted = sort(particles.col(d));
    
    for(int m = 0; m < M; m++){
      weights(m,d) = G(y(d), particles_sorted(m), sigma_y);
    }
    
    log_likelihood += log(mean(weights.col(d)));
    
    sum_weights = sum(weights.col(d));
    for(int m = 0; m < M; m++){
      weights(m,d) = weights(m,d)/sum_weights;
    }
  }
  
  return log_likelihood;
}