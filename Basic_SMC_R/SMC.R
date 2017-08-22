library(RSQMC) # Rcpp package
library(ggplot2)
set.seed(42)

# Parameters
T_ = 1000
sigma_y = 0.5**2
M = 10

# State-space model
# latent process is a random walk
x = matrix(0, T_)
x[1] = rnorm(1)
for(i in 2:T_){
  x[i] = rnorm(1, x[i-1])
}
y = rnorm(T_, x, sigma_y)

# Kernels and other functions
g <- function(y, x, sigma_2){
  return(dnorm(y, mean = x, sigma_2))
}

f <- function(x, x_prev){
  return(dnorm(x, mean = x_prev))
}

norm <- function(x){
  return (x/sum(x))
}

basic_SMC_R <- function(y, sigma_y, M, T_){
  # Filtering initialization
  particles = matrix(0, M, T_)
  W = matrix(0, M, T_)
  
  particles[, 1] = rnorm(M)
  
  W[,1] = g(y[1], particles[,1], sigma_y) * f(particles[, 1], 0)
  
  W[,1] = norm(W[,1])
  
  for(n in 2:T_){
    # Resampling
    resampling = sample(M, replace = T, prob = W[,n-1])
    
    # Mutation
    particles[, n] = rnorm(M, mean = particles[resampling, n-1])
    
    # Weights computation
    W[, n] = g(y[n], particles[, n], sigma_y)
    
    # Weights normalization
    W[, n] = norm(W[, n])
  }
  
  x_pred = matrix(0, T_)
  for(i in 1:T_){
    x_pred[i] = weighted.mean(particles[, i], W[, i])
  }
  
  # Backward smoothing/sampling
  backward_particles = matrix(0, M, T_)
  
  resampling = sample(M, replace = T, prob = W[, T_])
  
  backward_particles[, T_] = particles[resampling, T_]
  
  # Backward pass
  for(t in (T_-1):1){
    wgts = W[, t]*f(backward_particles[, t+1] ,particles[, t])
    resampling = sample(M, replace = T, prob = wgts)
    backward_particles[, t] = particles[resampling, t]
  }
  
  x_smooth = colMeans(backward_particles)
  return(x_smooth)
}

R_execution_time <- function(t){
  ptm <- proc.time()
  x_smooth <- basic_SMC_R(y, sigma_y, t, 30)
  return((proc.time() - ptm)[1])
}

Cpp_execution_time <-function(t){
  ptm_cpp <- proc.time()
  x_cpp = SMC(y, sigma_y, t, 30)
  return((proc.time() - ptm_cpp)[1])
}

R_elapsed = matrix(0, 500)
Cpp_elapsed = matrix(0, 500)
for(t in 20:500){
  
  R_elapsed[t] = mean(rep(R_execution_time(t), 50))
  
  # RCPP equivalent
  Cpp_elapsed[t] = mean(rep(Cpp_execution_time(t), 50))
}


df = data.frame(R=R_elapsed[20:500], Cpp=Cpp_elapsed[20:500])
df$time = 20:500

R_abline = coef(lm(R ~ time, data = df))
Cpp_abline = coef(lm(Cpp ~ time, data = df))

(ggplot(df) + geom_point(aes(x=time, y=R, color='R')) + 
  geom_point(aes(x=time, y=Cpp, color='C++')) + 
  geom_abline(intercept=R_abline[[1]], slope=R_abline[[2]], col='#00BFC4', size=1)+
  geom_abline(intercept=Cpp_abline[[1]], slope=Cpp_abline[[2]], col='red', size=1)+
  xlab('Number of particles') + ylab('Execution time') + ylim(0, 0.03))


#plot(x, col='blue')
#lines(x_cpp, col='black')
#lines(x_smooth, col='green')
