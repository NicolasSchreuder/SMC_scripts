library(RSQMC)
library(MASS)
library(mvtnorm)
source("/Users/Schreuder/Google Drive/ENSAE/2A/Stage2A/RSQMC_scripts/Lattice_Model_Functions.R") 

set.seed(42)

# SMC parameters
D <- 4 # dimension -> DxD square
N <- 100 # number of external particles
M <- 64 # number of internal particles
T_ = 15 # time

# AR(1) parameter
mu = 2
rho = 0.9

# Covariance matrix design
S_inv <- matrix(0, D**2, D**2)
tau = 1.
lambda = 1.

for(n in 1:D){
  for(k in 1:D){
    nb_neighbours = 0
    # (n,k) gives n*D + k position in vector
    vector_index <- (n-1)*D + k
    if (n > 1){ # upper neighbour
      S_inv[vector_index, vector_index - D] <- -lambda
      nb_neighbours = nb_neighbours + 1
    }
    if (n < D){ # lower neighbour
      S_inv[vector_index, vector_index + D] <- -lambda
      nb_neighbours = nb_neighbours + 1
    }
    if(k > 1){ # left neighbour
      S_inv[vector_index, vector_index - 1 ] <- -lambda 
      nb_neighbours = nb_neighbours + 1
    }
    if (k < D){ # right neighbour
      S_inv[vector_index, vector_index + 1] <- -lambda
      nb_neighbours = nb_neighbours + 1
    }
    S_inv[vector_index, vector_index] = tau + 2*nb_neighbours*lambda
  }
}

S <- solve(S_inv) # covariance matrix

# Calibration 
#S <- 0.5*S

# Trajectory simulation
x <- matrix(0, D**2, T_)
y <- matrix(0, D**2, T_)

x[, 1] <- mvrnorm(n = 1, matrix(mu, D**2), S)
for(d in 1:(D**2)){
  y[d, 1] <- rpois(1, exp(x[d, 1]))
}

for(t in 2:T_){
  x[, t] <- mu + rho*(x[, t-1]-mu) + mvrnorm(n = 1, matrix(0, D**2), S)
  for(d in 1:(D**2)){
    y[d, t] <- rpois(1, exp(x[d, t]))
  }
}

#### NSMC without QMC
particles <- lattice_filter(y, N, M, D, T_, mu, rho, S)

#### NSMC with QMC
particles_QMC <- lattice_filter_QMC(y, N, M, D, T_, mu, rho, S)

for(d in 1:D**2){
  plot(x[d, ])
  lines(x[d, ])
  lines(apply(particles[, d, 1:(T_-pred_horizon)], 2, mean), col='red')
  lines(apply(particles_QMC[, d, 1:(T_-pred_horizon)], 2, mean), col='blue')
}

#for(d in 1:D**2){
  #plot(x[d, ])
  #lines(x[d, ])
  #polygon(c(1:(T_-1), rev(1:(T_-1))), c(apply(particles[, d, 1:(T_-pred_horizon)], 2, quantile)[2, ], rev(apply(particles[, d, 1:(T_-pred_horizon)], 2, quantile)[4, ])), col=rgb(1,0, 1,0.1), border = NA)
  #lines(apply(particles[, d, 1:(T_-pred_horizon)], 2, quantile)[3, ], col='blue')
#}

### Prediction

pred_particles = array(0, dim=c(N, D**2, pred_horizon))
y_pred = array(0, dim=c(D**2, pred_horizon))

for(n in 1:N){
  pred_particles[n, , 1] = mu + rho*(particles[n, , T_-1] - mu) + mvrnorm(n = 1, matrix(0, D**2), S)
}

x_quantile = apply(pred_particles[ , , 1], 2, quantile)

for(d in 1:D**2){
  plot(x[d, ])
  lines(x[d, ])
  lines(apply(particles[, d, 1:(T_-pred_horizon)], 2, mean), col='red')
  points(T_, x_quantile[2, d], col='green')
  points(T_, x_quantile[4, d], col='green')
  points(T_, x_quantile[3, d], col='blue')
}