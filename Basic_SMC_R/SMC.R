set.seed(42)

# Parameters
T = 20
sigma = 0.5
M = 200

# State-space model
x = matrix(0, T)
x[1] = rnorm(1)
for(i in 2:T){
  x[i] = rnorm(1,x[i-1])
}
y = rnorm(T,x,sigma)

# Kernels and other functions
g <- function(y,x){
  return(dnorm(y, mean = x, sigma))
}

f <- function(x,x_prev){
  return(dnorm(x, mean = x_prev))
}

norm <- function(x){
  return (x/sum(x))
}

ptm <- proc.time()
# Filtering initialization
particles = matrix(0, M, T)
W = matrix(0, M, T)

particles[,1] = rnorm(M)

W[,1] = g(y[1], particles[,1]) * f(particles[,1], 0)

W[,1] = norm(W[,1])

for(n in 2:T){
  # Resampling
  resampling = sample(M, replace = T, prob = W[,n-1])

  # Mutation
  particles[,n] = rnorm(M, mean = particles[resampling, n-1])
  
  # Weights computation
  W[,n] = g(y[n], particles[,n])
  
  # Weights normalization
  W[,n] = norm(W[,n])
}

x_pred = matrix(0,T)
for(i in 1:T){
  x_pred[i] = weighted.mean(particles[,i], W[,i])
}

# Backward smoothing/sampling
backward_particles = matrix(0, M, T)

resampling = sample(M, replace = T, prob = W[,T])

backward_particles[,T] = particles[resampling,T]

# Backward pass
for(t in (T-1):1){
    wgts = W[,t]*f(backward_particles[,t+1] ,particles[,t])
    resampling = sample(M, replace = T, prob = wgts)
    backward_particles[,t] = particles[resampling, t]
}

x_smooth = colMeans(backward_particles)

print(proc.time() - ptm)

# RCPP equivalent
ptm_cpp <- proc.time()
x_cpp = SMC(y, M, T)
print(proc.time() - ptm_cpp)

plot(x, col='blue')
lines(x)
lines(x_pred, col='red')

#lines(x_cpp, col='black')
#lines(x_smooth, col='green')


