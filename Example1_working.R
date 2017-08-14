library(MASS) # for multivariate Gaussian simulation
library(RSQMC) # for internal filters
library(profvis) # R profiler
source("/Users/Schreuder/Google Drive/ENSAE/2A/Stage2A/bRistol/Example1Functions.R") 

set.seed(42)

## Parameters ##
T_ <- 40 # time
D <- 64 # dimension of x
sigma_y <- .25**2 # variance of y|x
M <- 16 # number of internal particles
N <- 100 # number of external particles
tau <- 1. # coordinates variance
lambda <- 1. # coordinates correlation
nb_filters <- 1 # number of external filters 

## Covariance matrix 
S = covariance_matrix(D, tau, lambda)

# State-space model
SSM = state_space_model(D, T_, S, sigma_y)
x = SSM$X
y = SSM$Y

## Without QMC ## 

profvis({
pred_without_QMC <- array(0, dim=c(nb_filters, D, T_))

particles <- array(0, dim=c(N, D, T_+1))
weights <- matrix(0, N)
log_likelihood <- matrix(0, N)
internal_weights <- array(0, dim=c(T_+1, M, D, N))
internal_particles <- array(0, dim=c(T_+1, M, D, N))

for(f in 1:nb_filters){
  print(c("Filter without QMC n°", f))

  # External particles initialization
  particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)
  
  for(t in 2:(T_+1)){
    # Internal filter 
    # Returns the log-likelihood ( = weights for resampling)
    internal_filter_output = internal_filter(N, y[, t-1], M, D, S, particles[, , t-1], sigma_y)
    internal_weights[t, , , ] <- internal_filter_output$Weights
    internal_particles[t, , ,] <- internal_filter_output$Particles
    log_likelihood <- internal_filter_output$Log_Likelihood
  
    
    # Weights normalization
    weights = exp(log_likelihood - max(log_likelihood))
    sum = sum(weights)
    weights = weights/sum

    # Resampling
    resampling <- sample(N, replace = T, prob = weights)
    
    # Backward_particles
    sampling_weights <- matrix(0,M)
    
    for(n in 1:N){
      component_resampling <- sample(M, 1, replace = T, internal_weights[t, , D, resampling[n]])
      particles[n, D, t] <-  internal_particles[t, component_resampling, D, resampling[n]]
      for(d in D:2){
        sd <- S[d, d] - (S[d-1, d]**2)/S[d-1, d-1];
        mu <- 0.5*particles[resampling[n], d, t-1] + S[d-1, d]/S[d-1, d-1]*(internal_particles[t, , d-1, resampling[n]] - 0.5*particles[resampling[n], d-1, t-1])
        sampling_weights <- internal_weights[t, , d-1, resampling[n]] * dnorm(internal_particles[t, , d, resampling[n]], mu, sd)
        component_resampling <- sample(M, 1, replace = T, sampling_weights)
        particles[n, d-1, t] <- internal_particles[t, component_resampling, d-1, resampling[n]]
      }
    }
  }
  pred_without_QMC[f, , ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
}
})



### QMC ###

pred_with_QMC <- array(0, dim=c(nb_filters, D, T_))

particles <- array(0, dim=c(N, D, T_+1))
weights <- matrix(0, N)
log_likelihood <- matrix(0, N)
internal_weights <- array(0, dim=c(T_+1, M, D, N))
internal_particles <- array(0, dim=c(T_+1, M, D, N))

for(f in 1:nb_filters){
  print(c("Filter with QMC", f))
  
  # Initialize particles : particles * components * time
  particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)
  
  for(t in 2:(T_+1)){
    internal_filter_output = internal_filter_with_QMC(N, y[, t-1], M, D, S, particles[, , t-1], sigma_y)
    internal_weights[t, , , ] <- internal_filter_output$Weights
    internal_particles[t, , ,] <- internal_filter_output$Particles
    log_likelihood <- internal_filter_output$Log_Likelihood
    
    # Weights normalization
    weights = exp(log_likelihood - max(log_likelihood))
    sum = sum(weights)
    weights = weights/sum
    
    # Resampling
    resampling = sample(N, replace = T, prob = weights)
  
    sampling_weights <- matrix(0,M)
    
    for(n in 1:N){
      component_resampling <- sample(M, 1, replace = T, internal_weights[t, , D, resampling[n]])
      particles[n, D, t] <-  internal_particles[t, component_resampling, D, resampling[n]]
      for(d in D:2){
        sd <- S[d, d] - (S[d-1, d]**2)/S[d-1, d-1];
        mu <- 0.5*particles[resampling[n], d, t-1] + S[d-1, d]/S[d-1, d-1]*(internal_particles[t, , d-1, resampling[n]] - 0.5*particles[resampling[n], d-1, t-1])
        sampling_weights <- internal_weights[t, , d-1, resampling[n]] * dnorm(internal_particles[t, , d, resampling[n]], mu, sd)
        component_resampling <- sample(M, 1, replace = T, sampling_weights)
        particles[n, d-1, t] <- internal_particles[t, component_resampling, d-1, resampling[n]]
      }
    }
  }
  pred_with_QMC[f, , ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
}

library(ggplot2)
plot_coordinate = 2
df = data.frame(real_value = x[plot_coordinate, ], no_QMC=pred_without_QMC[1, plot_coordinate, ], QMC=pred_with_QMC[1, plot_coordinate, ])
df$idu <- as.numeric(row.names(df))
ggplot(df) + geom_point(aes(x=idu, y=real_value, color='Real value'))+ geom_line(aes(x=idu, y=no_QMC, color='Without QMC')) + geom_line(aes(x=idu, y=QMC, color='With QMC')) 

var = apply(pred_without_QMC[ , , ], c(2,3), var)
mean_var= apply(var , 2, mean)

var_QMC = apply(pred_with_QMC[ , ,], c(2,3), var)
mean_var_QMC = apply(var_QMC , 2, mean)

var_ratio = data.frame(variance_ratio = apply(var_QMC/mean_var_QMC, 2, mean))
ggplot(var_ratio) + geom_line(aes(x=1:T_, y=variance_ratio), colour='blue') + geom_hline(yintercept = 1)

if (0==1){
## Naive SMC ##
library(mvtnorm)

pred_naive <- array(0, dim=c(D, T_))

# External particles initialization
particles <- array(0, dim=c(N, D, T_+1))
particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)

W <- matrix(0, N, T_)

for(n in 1:N){
  W[n, 1] = dmvnorm(x=y[, 1], mean=particles[n, ,1], sigma=diag(sigma_y, D)) * dmvnorm(x=particles[n, ,1], mean=matrix(0, D), sigma=S)
}

W[, 1] = W[,1]/sum(W[, 1])

for(t in 2:(T_)){
  # Resampling
  resampling <- sample(N, replace = T, prob = W[, t-1])
  
  for(n in 1:N){
    particles[n, ,t] = mvrnorm(n = 1, particles[resampling[n] , ,t-1], S)
  }
  
  for(n in 1:N){
    W[n, t] = dmvnorm(x=y[, t], mean=particles[n, , t], sigma=diag(sigma_y, D))
  }
  
  W[, t] = W[,t]/sum(W[, t])
  
}
pred_naive[, ] <- apply(particles[ , ,2:(T_)], 2, mean)

plot_coordinate = 1
plot(x[plot_coordinate, ])
lines(x[plot_coordinate, ], col="black")
lines(pred_naive[plot_coordinate, ], col="green")
}