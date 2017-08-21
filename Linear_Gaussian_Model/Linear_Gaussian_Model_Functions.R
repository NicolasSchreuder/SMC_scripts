library(MASS) # for multivariate Gaussian simulation

covariance_matrix <- function(D, tau, lambda){
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
  return(S)
}

# Generates latent variable and observation variable from the SSM
state_space_model <- function(D, T_, S, sigma_y){
  x <- matrix(0, D, T_)
  y <- matrix(0, D, T_)
  
  x[, 1] <- mvrnorm(n = 1, matrix(0, D), S)
  y[, 1] <- mvrnorm(n = 1, x[, 1], diag(sigma_y, D))
  
  for(t in 2:T_){
    x[, t] <- mvrnorm(n = 1, 0.5 * x[, t-1], S)
    y[, t] <- mvrnorm(n = 1, x[, t],  diag(sigma_y, D))
  }
  return(list(X=x,Y=y))
}

backward_paralell <- function(n, resampling, N, M, t, D, weights, particles, internal_particles, internal_weights, S){
  particle = matrix(0, D)
  sampling_weights <- matrix(0,M)
  
  component_resampling <- sample(M, 1, replace = T, internal_weights[t, , D, resampling[n]])
  particles[D] <-  internal_particles[t, component_resampling, D, resampling[n]]
  
  for(d in D:2){
    sd <- S[d, d] - (S[d-1, d]**2)/S[d-1, d-1];
    mu <- 0.5*particles[resampling[n], d, t-1] + S[d-1, d]/S[d-1, d-1]*(internal_particles[t, , d-1, resampling[n]] - 0.5*particles[resampling[n], d-1, t-1])
    sampling_weights <- internal_weights[t, , d-1, resampling[n]] * dnorm(internal_particles[t, , d, resampling[n]], mu, sd)
    component_resampling <- sample(M, 1, replace = T, sampling_weights)
    particle[d-1] <- internal_particles[t, component_resampling, d-1, resampling[n]]
  }
  return(particle)
}

# Bootstrap filter for the gaussian linear model
gaussian_bootstrap_filter <- function(y, sigma_y, D, T_, S, number_of_particles){
  
  # Particles initialization
  particles <- array(0, dim=c(number_of_particles, D, T_+1))
  particles[, , 1] <- mvrnorm(n = number_of_particles, matrix(0, D), S)
  
  # Weights initialization
  W <- matrix(0, number_of_particles, T_)
  
  # Initial importance weights
  for(n in 1:number_of_particles){
    W[n, 1] = dmvnorm(y[, 1], mean=particles[n, , 1], sigma=diag(sigma_y, D))
  }
  
  # Weights normalization
  W[, 1] = W[, 1]/sum(W[, 1])
  
  for(t in 2:T_){
    # Resampling
    resampling <- sample(number_of_particles, replace = T, prob = W[, t-1])
    
    # Propagation
    for(n in 1:number_of_particles){
      particles[n, ,t] = mvrnorm(n = 1, particles[resampling[n], , t-1], S)
    }
    
    # Importance weights
    for(n in 1:number_of_particles){
      W[n, t] = dmvnorm(y[, t], mean=particles[n, , t], sigma=diag(sigma_y, D))
    }
    
    # Weights normalization
    W[, t] = W[, t]/sum(W[, t])
  }
  
  # Backward smoothing/sampling
  backward_particles = array(0, dim=c(number_of_particles, D, T_))
  
  resampling = sample(number_of_particles, replace = T, prob = W[, T_])
  backward_particles[, , T_] = particles[resampling, , T_]
  
  backward_weights = array(0, N)
  
  # Backward pass
  for(t in (T_-1):1){
    
    # Backward weights
    for(n in 1:number_of_particles){
      backward_weights[n] = W[n, t] * dmvnorm(backward_particles[n, , t+1], particles[n, , t], S)
    }
    
    resampling = sample(number_of_particles, replace = T, prob = backward_weights)
    backward_particles[, , t] = particles[resampling, , t]
  }
  
  # Latent variable MAP
  pred_naive <- array(0, dim=c(D, T_))
  
  pred_naive[, 2:T_] <- apply(backward_particles[, , 2:(T_)], c(2, 3), mean)
  
  return(pred_naive)
}


gaussian_guided_filter <- function(y, sigma_y, D, T_, S, number_of_particles){
  # Optimal proposal covariance matrix
  A <- solve(S) + diag(1/sigma_y, D)
  
  # External particles initialization
  particles <- array(0, dim=c(number_of_particles, D, T_+1))
  
  b <- 0.5*solve(A)%*%((2/sigma_y)*y[, 1])
  particles[, , 1] <- mvrnorm(n = number_of_particles, b, solve(A))
  
  # Resampling weights are all equal to 1/N when using the optimal proposal
  resampling_weights = matrix(1/number_of_particles, number_of_particles)
  
  for(t in 2:T_){
    # Resampling
    resampling <- sample(number_of_particles, replace = T, prob = resampling_weights)
    
    # Propagation
    for(n in 1:number_of_particles){
      # Optimal proposal mean
      b <- 0.5*solve(A)%*%(solve(S)%*%particles[resampling[n], , t-1] + (2/sigma_y)*y[, t]) 
      particles[n, , t] = mvrnorm(n = 1, b, solve(A))
    }
  }
  
  # Backward smoothing/sampling
  backward_particles = array(0, dim=c(number_of_particles, D, T_))
  
  resampling = sample(number_of_particles, replace = T_, prob = resampling_weights)
  
  backward_particles[, , T_] = particles[resampling, , T_]
  
  backward_weights = array(0, number_of_particles)
  
  # Backward pass
  for(t in (T_-1):1){
    for(n in 1:number_of_particles){
      backward_weights[n] = 1/number_of_particles*dmvnorm(backward_particles[n, , t+1], particles[n, , t], S)
    }
    resampling = sample(number_of_particles, replace = T, prob = backward_weights)
    backward_particles[, , t] = particles[resampling, , t]
  }
  
  pred <- array(0, dim=c(D, T_))
  pred[, 2:T_] <- apply(particles[, , 2:T_], c(2,3) , mean)
  return(pred)
}



# NSMC for the linear gaussian model
NSMC_linear_gaussian <- function(y, N, M, D, T_, S, sigma_y, QMC){
  
  particles <- array(0, dim=c(N, D, T_+1))
  weights <- array(0, N)
  log_likelihood <- array(0, N)
  internal_weights <- array(0, dim=c(T_+1, M, D, N))
  internal_particles <- array(0, dim=c(T_+1, M, D, N))
  
  # External particles initialisation 
  particles[, , 1] <- mvrnorm(n = N, matrix(0, D), S)
  
  for(t in 2:(T_+1)){
    
    # The internal filter returns :
    # - the log-likelihood of the SMC along the coordinates (= weights for resampling) 
    # - the generated internal particles (needed for backward simulation)
    # - the associated internal weights (needed for backward simulation)
    # NB : the particles start at time t=2 and the observations y at t=1 
    if(QMC){
      internal_filter_output = internal_filter_with_QMC(N, y[, t-1], M, D, S, 
                                                        particles[, , t-1], sigma_y)
    } else {
      internal_filter_output = internal_filter(N, y[, t-1], M, D, S, 
                                               particles[, , t-1], sigma_y)
    }
    
    internal_weights[t, , , ] <- internal_filter_output$Weights
    internal_particles[t, , ,] <- internal_filter_output$Particles
    log_likelihood <- internal_filter_output$Log_Likelihood
    
    # Weights normalisation
    weights = exp(log_likelihood - max(log_likelihood))
    sum = sum(weights)
    weights = weights/sum
    
    # Resampling step
    resampling = sample(N, replace = T, prob = weights)
    
    # Backward simulation (parallelised)
    cl <- makeCluster(3)
    particles[, , t] <- t(parSapply(cl, 1:N, backward_paralell, resampling=resampling, 
                                    N=N, M=M, t=t, D=D, weights=weights, particles=particles, 
                                    internal_particles=internal_particles, 
                                    internal_weights=internal_weights, S=S))
    stopCluster(cl)
  }
  return(particles[ , , 2:(T_+1)])
}