# Filtre externe temporel + filtre interne sur les coordonnées
library(MASS) # for multinormal simulation
library(RSQMC)
set.seed(42)

# Parameters
T_ <- 50 # time
D <- 16 # dimension of x
sigma_y <- .25**2 # variance of y|x
M <- 512 # number of internal particles
N <- 200 # number of external particles
tau <- 1. # coordinates variance
lambda <- 1. # coordinates correlation
nb_filters <- 1 # number of external filters 

# Covariance matrix
S_inv <- matrix(0, D, D)
S_inv[1, 1] <- tau + lambda
S_inv[D, D] <- tau + lambda

for(i in 2:(D-1)){
  S_inv[i, i] <- tau + 2*lambda
}

for(i in 1:(D-1)){
  S_inv[i, i+1] <- -lambda
  S_inv[i+1, i] <- -lambda
}

S <- solve(S_inv)

# State-space model
x <- matrix(0, D, T_)
y <- matrix(0, D, T_)

x[, 1] <- mvrnorm(n = 1, matrix(0, D), S)
y[, 1] <- mvrnorm(n = 1, x[, 1], diag(sigma_y, D))

for(t in 2:T_){
  x[, t] <- mvrnorm(n = 1, 0.5 * x[, t-1], S)
  y[, t] <- mvrnorm(n = 1, x[, t],  diag(sigma_y, D))
}

library(profvis)
profvis({
pred_without_QMC <- array(0, dim=c(nb_filters, D, T_))
  
particles <- array(0, dim=c(N, D, T_+1))
weights <- matrix(0, N)
log_likelihood <- matrix(0, N)
internal_weights <- array(0, dim=c(N, T_+1, M, D))
internal_particles <- array(0, dim=c(N, T_+1, M, D))

for(f in 1:nb_filters){
  print(c("Filter without QMC n°", f))
  
  # External particles initialization
  
  particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)

  for(t in 2:(T_+1)){
    # internal filter for each external particle 
    # returns the log-likelihood ( = weights for resampling)
    for(n in 1:N){
      internal_filter_output <- SMC_example_1(y[,t-1], M, D, S, particles[n, , t-1], sigma_y)
      internal_weights[n, t, , ] <- internal_filter_output$Weights
      internal_particles[n, t, ,] <- internal_filter_output$Particles
      log_likelihood[n] <- internal_filter_output$Log_Likelihood
    }
    
    # Weights normalization
    weights = exp(log_likelihood - max(log_likelihood))
    sum = sum(weights)
    weights = weights/sum

    # Resampling
    resampling <- sample(N, replace = T, prob = weights)
    
    # Mutation
    #variance <- solve(S_inv + diag(1/sigma_y, D))
    #for(n in 1:N){
    #  mu <- 0.5 * variance %*% (S_inv %*% particles[resampling[n], , t-1] + (2/sigma_y)*y[, t-1])
    #  particles[n, , t] <- mvrnorm(n = 1, mu, variance)
    #}
    
    # Backward_particles
    for(n in 1:N){
      component_resampling <- sample(M, 1, replace = T, internal_weights[resampling[n], t, , D])
      particles[n, D, t] <-  internal_particles[resampling[n], t, component_resampling, D]
      for(d in D:2){
        sd <- S[d, d] - (S[d-1, d]**2)/S[d-1, d-1];
        sampling_weights <- matrix(0,M)
        mu <- 0.5*particles[resampling[n], d, t-1] + S[d-1, d]/S[d-1, d-1]*(internal_particles[resampling[n], t, , d-1] - 0.5*particles[resampling[n], d-1, t-1])
        sampling_weights <- internal_weights[resampling[n], t, , d-1] * dnorm(particles[resampling[n], d, t], mu, sd)
        component_resampling <- sample(M, 1, replace = T, sampling_weights)
        particles[n, d-1, t] <- internal_particles[resampling[n], t, component_resampling,d-1]
      }
    }
  }
  pred_without_QMC[f, , ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
}
})

pred_without_QMC <- array(0, dim=c(D, T_+1))
pred_without_QMC[, ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
plot(x[3,])
lines(x[3,])
lines(pred_without_QMC[1, 3, ], col="green")

### QMC ###

pred_with_QMC = matrix(0, nb_filters, T_+1)

for(f in 1:nb_filters){
  print(c("Filter with QMC", f))
  
  # Initialize particles : particles * components * time
  particles <- array(0, dim=c(N, D, T_+1))
  particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)
  
  weights = matrix(0, N)
  log_likelihood = matrix(0, N)

  for(t in 2:(T_+1)){

    # Internal filter for each particle -> returns the log-likelihood
    for(n in 1:N){
      internal_filter_output <- SMC_example_1_QMC(y[,t-1], M, D, S, particles[n, , t-1], sigma_y)
      internal_weights[n, t, , ] <- internal_filter_output$Weights
      internal_particles[n, t, ,] <- internal_filter_output$Particles
      log_likelihood[n] <- internal_filter_output$Log_Likelihood
    }
    
    # Weights normalization
    weights = exp(log_likelihood - max(log_likelihood))
    sum = sum(weights)
    weights = weights/sum
    
    # Resampling
    resampling = sample(N, replace = T, prob = weights)
    
    # Mutation
    #variance = solve(S_inv + diag(1/sigma_y, D))
    #for(n in 1:N){
    #  mu = 0.5 * variance %*% (S_inv %*% particles[resampling[n], , t-1] + (2/sigma_y)*y[,t-1])
    #  particles[n, , t] = mvrnorm(n=1, mu, variance)
    #}
    
    # Backward_particles
    for(n in 1:N){
      component_resampling <- sample(M, 1, replace = T, internal_weights[resampling[n], t, , D])
      particles[n, D, t] <-  internal_particles[resampling[n], t, component_resampling, D]
      for(d in D:2){
        sd <- S[d, d] - (S[d-1, d]**2)/S[d-1, d-1];
        sampling_weights <- matrix(0,M)
        mu <- 0.5*particles[resampling[n], d, t-1] + S[d-1, d]/S[d-1, d-1]*(internal_particles[resampling[n], t, , d-1] - 0.5*particles[resampling[n], d-1, t-1])
        sampling_weights <- internal_weights[resampling[n], t, , d-1] * dnorm(particles[resampling[n], d, t], mu, sd)
        
        component_resampling <- sample(M, 1, replace = T, sampling_weights)
        particles[n, d-1, t] <- internal_particles[resampling[n], t, component_resampling,d-1]
      }
    }
  }
}

pred_with_QMC <- array(0, dim=c(D, T_+1))
pred_with_QMC[, ] <- apply(particles[ , ,], c(2,3), mean)
lines(pred_with_QMC[3, ], col="green")

var_err = apply(pred_without_QMC, 2, var)
var_err_QMC = apply(pred_with_QMC, 2, var)
plot(var_err/var_err_QMC, type='l')
abline(a=1, b=0, col='red')

x = gen_sobol(2, 256)
# plot(x[1,], x[2, ])