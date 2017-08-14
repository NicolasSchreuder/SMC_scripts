# Example #
library(MASS) # for multinormal simulation

# Parameters
T = 50
sigma = 0.5
M = 50

# Covariance matrix
tau = 0.1
lambda = 0.5
S_inv = matrix(0,T,T)
S_inv[1,1] = tau + lambda
S_inv[T,T] = tau + lambda
for(i in 2:(T-1)){
  S_inv[i,i] = tau + 2*lambda
}

for(i in 1:T-1){
  S_inv[i, i+1] = -lambda 
  S_inv[i+1, i] = -lambda
}

S = solve(S_inv)

# State-space model
# Equivalent de x(t-1)
zeros = matrix(0, T)
x_init = mvrnorm(n=1, zeros, S)

x = matrix(0, T)
x[1] = rnorm(1)
for(i in 2:T){
  x[i] = rnorm(1, 0.5*x_init[i] + S[i-1,i]/S[i,i]*(x[i-1] - 0.5*x_init[i-1]), S[i,i] - (S[i-1,i]**2)/S[i-1,i-1])
}
y = rnorm(T,x,sigma)

# Does not change 
g <- function(y,x){
  return(dnorm(y, mean = x, sigma))
}

# x(k) dépend de x(t-1, k) (donné) et de x(t,k-1)
f <- function(x,x_prev, x_init, i){
  return(dnorm(x, mean = 0.5*x_init[i] + S[i-1,i]/S[i,i]*(x[i-1] - 0.5*x_init[i-1]), S[i,i] - (S[i-1,i]**2)/S[i-1,i-1] ))
}

norm <- function(x){
  return (x/sum(x))
}

ptm <- proc.time()
# Filtering initialization
particles = matrix(0, M, T)
W = matrix(0, M, T)

particles[,1] = rnorm(M, mean=0, S[1,1])

W[,1] = g(y[1], particles[,1]) * dnorm(particles[,1], mean = 0.5*x_init[1], S[1,1] )

W[,1] = norm(W[,1])

for(n in 2:T){
  # Resampling
  resampling = sample(M, replace = T, prob = W[,n-1])
  
  mu = 0.5*x_init[n] + S[n-1,n]/S[n,n]*(particles[resampling, n-1] - 0.5*x_init[n-1]) 
  
  sd = S[n,n] - (S[n-1,n]**2)/S[n-1,n-1]
    
  # Mutation
  particles[,n] = rnorm(M, mu, sd)
  
  # Weights computation
  W[,n] = g(y[n], particles[,n])
  
  # Weights normalization
  W[,n] = norm(W[,n])
}

x_pred = matrix(0,T)
for(i in 1:T){
  x_pred[i] = weighted.mean(particles[,i], W[,i])
}
print(proc.time() - ptm)
# Backward smoothing/sampling
backward_particles = matrix(0, M, T)

resampling = sample(M, replace = T, prob = W[,T])

backward_particles[,T] = particles[resampling,T]

# Backward pass
for(t in (T-1):1){
  # Il faut revoir les arguments de f ici
  wgts = W[,t]*f(backward_particles[,t+1], 0, particles[,t], t+1)
  resampling = sample(M, replace = T, prob = wgts)
  backward_particles[,t] = particles[resampling, t]
}

x_smooth = colMeans(backward_particles)


# RCPP equivalent
ptm_cpp <- proc.time()
x_cpp = SMC_example_1(y, M, T, S, x_init)
print(proc.time() - ptm_cpp)


plot(x, col='blue')
lines(x_pred, col='red')
lines(x_cpp, col='green')
lines(x_smooth, col='black')

