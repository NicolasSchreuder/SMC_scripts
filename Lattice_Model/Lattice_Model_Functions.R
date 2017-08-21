# Multivariate density of independant Poisson variables
observation_lh <- function(y, x){
  prod <- 1
  for(i in 1:length(y)){
    prod <- prod * dpois(y[i], exp(x[i]))
  }
  return(prod)
}

lattice_filter <- function(y, N, M, D, T_, mu, rho, S, QMC=FALSE){
  
  particles <- array(0, dim=c(N, D**2, T_))
  resampling_weights <- array(0, N)
  
  for(n in 1:N){
    particles[n, , 1] <- mvrnorm(n = 1, matrix(mu, D**2), S)
    resampling_weights[n] <- observation_lh(y[, 1], particles[n, , 1])
  }
  
  resampling_weights = resampling_weights/sum(resampling_weights)
  
  resampling = sample(N, N, replace=T, resampling_weights)
  
  a <- particles[, , 1]
  particles[, , 1] <- a[resampling, ]
  
  Z = matrix(1, N)
  internal_particles <- array(0, dim=c(M, D**2, N))
  
  for(t in 2:T_){
    
    if(QMC){
      internal_filter_output = lattice_internal_filter_QMC(N, M, D, 
                                                       mu + rho*(particles[, , t-1] - mu), 
                                                       S, y[, t])
    } else {
      internal_filter_output = lattice_internal_filter(N, M, D, 
                                                           mu + rho*(particles[, , t-1] - mu), 
                                                           S, y[, t])
    }
    
    log_lh = internal_filter_output$Log_likelihood
    Z = exp(log_lh - min(log_lh))
    Z = Z/sum(Z)
    
    internal_particles = internal_filter_output$Particles
    
    resampling_ancestor = sample(N, N, replace=T, Z)
    sampling = sample(M, N, replace=T, matrix(1/M, M))
    
    for(n in 1:N){
      particles[n, , t] <- internal_particles[sampling[n], , resampling_ancestor[n]]
    }
  }
  return(particles)
}

lattice_filter_replicate <-function(i){
  particles <- array(0, dim=c(N, D**2, T_))
  particles <- lattice_filter(y, N, M, D, T_, mu, rho, S, QMC=FALSE)
  return(apply(particles[, 1, ], 2, mean))
}

lattice_filter_QMC_replicate <-function(i){
  particles <- array(0, dim=c(N, D**2, T_))
  particles <- lattice_filter(y, N, M, D, T_, mu, rho, S, QMC=TRUE)
  return(apply(particles[, 1, ], 2, mean))
}