library(MASS) # for multivariate Gaussian simulation
library(RSQMC) # for internal filters
library(profvis) # R profiler
library(parallel)
library(ggplot2)
source("/Users/Schreuder/Google Drive/ENSAE/2A/Stage2A/bRistol/Example1Functions.R") 

set.seed(41)

## Parameters ##
T_ <- 20 # time
D <- 8 # dimension of x
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

df = data.frame(x = x[2, ], y= y[2,])
df$idu <- as.numeric(row.names(df))
ggplot(df) + geom_point(aes(x=idu, y=y, color='Observation'), size=2)+ geom_line(aes(x=idu, y=x, color='Latent variable'), size=1) + xlab('Time') + ylab('Value')

## Without QMC ## 

#profvis({
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
    print(t)
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

    resampling = sample(N, replace = T, prob = weights)
    
    # Backward
    cl <- makeCluster(3)
    particles[, , t] <- t(parSapply(cl, 1:N, backward_paralell, resampling=resampling, N=N, M=M, t=t, D=D, weights=weights, particles=particles, internal_particles=internal_particles, internal_weights=internal_weights, S=S))
    stopCluster(cl)
  }
  pred_without_QMC[f, , ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
}
#})

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
    print(t)
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
    
    cl <- makeCluster(3)
    particles[, , t] <- t(parSapply(cl, 1:N, backward_paralell, resampling=resampling, N=N, M=M, t=t, D=D, weights=weights, particles=particles, internal_particles=internal_particles, internal_weights=internal_weights, S=S))
    stopCluster(cl)
  }
  pred_with_QMC[f, , ] <- apply(particles[ , ,2:(T_+1)], c(2,3), mean)
}


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

vars = data.frame(var = apply(var, 2, mean), var_QMC = apply(var_QMC, 2, mean))
ggplot(vars) + geom_line(aes(x=1:T_, y=var, colour='blue'), size=1.5) + geom_line(aes(x=1:T_, y=var_QMC, colour='green'), size=1.5) + scale_color_discrete(name = "Variances", labels = c("Without QMC", "With QMC"))


if (1==1){
## Naive SMC ## Without optimal proposal
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

# Backward smoothing/sampling
backward_particles = array(0, dim=c(N, D, T_))

resampling = sample(N, replace = T_, prob = W[,T_])

backward_particles[,,T_] = particles[resampling, ,T_]

wgts = array(0, N)

# Backward pass
for(t in (T_-1):1){
  for(n in 1:N){
  wgts[n] = W[n,t]*dmvnorm(backward_particles[n, ,t+1], particles[n, , t], S)
  }
  resampling = sample(N, replace = T, prob = wgts)
  backward_particles[,,t] = particles[resampling, , t]
}

pred_naive[2, 2:T_] <- apply(backward_particles[, 2,2:(T_)], 2, mean)

plot_coordinate = 2

naive = data.frame(real_value = x[plot_coordinate, ], pred = pred_naive[plot_coordinate, ])
ggplot(naive) + geom_line(aes(x=1:T_, y=real_value, colour='Real value'), size=1) + geom_line(aes(x=1:T_, y=pred, colour='Prediction'), size=1) + xlab('Time') + ylab('Value') + annotate("text", x = 18, y = 1, label = "n_x = 2")
}