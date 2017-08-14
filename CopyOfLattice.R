library(mvtnorm)
library(RSQMC)
source("/Users/Schreuder/Google Drive/ENSAE/2A/Stage2A/bRistol/Example1Functions.R") 

# SMC parameters
D <- 2 # dimension -> DxD square
N <- 50 # number of external particles
M <- 16 # number of internal particles
T_ = 15 # time
sigma_y = 0.25**2

mu = 1
rho = 0.9

# Multivariate density of independant Poisson variables
g <- function(y, x){
  prod <- 1
  for(i in 1:length(y)){
    prod <- prod * dpois(y[i], exp(x[i]))
    #prod <- prod * dnorm(y[i], x[i], sigma_y)
  }
  return(prod)
}


# Covariance matrix construction
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

S <- 4*S

# Trajectory simulation
x <- matrix(0, D**2, T_)
y <- matrix(0, D**2, T_)

x[, 1] <- mvrnorm(n = 1, matrix(mu, D**2), S)
for(d in 1:(D**2)){
  y[d, 1] <- rpois(1, exp(x[d, 1]))
  #y[d, 1] <- rnorm(1, x[d, 1], sigma_y)
}

for(t in 2:T_){
  x[, t] <- mu + rho*(x[, t-1]-mu) + mvrnorm(n = 1, matrix(0, D**2), S)
  for(d in 1:(D**2)){
    y[d, t] <- rpois(1, exp(x[d, t]))
    #y[d, t] <- rnorm(1, x[d, t], sigma_y)
  }
}

# 2^{nb_etapes} = D^2 =>
nb_etapes <- log(D**2)/log(2)

particles <- array(0, dim=c(N, D**2, T_))
resampling_weights <- array(0, N)

# Peut etre qu'il faut lancer un filtre interne pour l'initialisation ?
for(p in 1:N){
  particles[p, , 1] <- mvrnorm(n = 1, matrix(mu, D**2), S)
  resampling_weights[p] <- g(y[, 1], particles[p, , 1])
}

resampling_weights = resampling_weights/sum(resampling_weights)

resampling = sample(N, N, replace=T, resampling_weights)

a <- particles[, , 1]
particles[, , 1] <- a[resampling, ]

for(t in 2:T_){
  print(t)
  # Filtre interne lancé pour chaque particule externe
  internal_particles <- array(0, dim=c(N, M, D**2))
  updated_internal_particles <- array(0, dim=c(N, M, D**2))
  weights <- array(0, dim=c(M, D**2))
  
  # Likelihood estimate
  Z = matrix(1, N)
  
  if(F){
    for(p in 1:N){
      x_prev <- mu + rho*(particles[p, , t-1] - mu)
      
      for(k in 0:nb_etapes){
        print(k)
        d = 1 # indice de composante
        while(d <= D**2){
  
          if(k==0){ # Etape 0 : initalisation
            for(m in 1:M){
              # On initialise chaque particule en simulant selon la loi marginale
              internal_particles[p, m, d] <- x_prev[d] + rnorm(1, 0, S[d, d])
              weights[m, d] <- g(y[d, t], internal_particles[p, m, d])
            }
          }
          
          if(k==1){ # Etape 1 : on regroupe par couple (x_k, x_{k+1})
            # on regarde toutes les combinaisons possibles de particules et on cacule les poids associés
            v <- array(1, dim=c(M, M))
            for(n in 1:M){
              for(m in 1:M){
                v[n, m] <- ((weights[n, d]*weights[m, d+1])*  
                  dmvnorm(c(internal_particles[p, n, d], internal_particles[p, m, d+1]), x_prev[d:(d+1)], S[d:(d+1), d:(d+1)]) /
                  (dnorm(internal_particles[p, n, d], x_prev[d], S[d, d])*dnorm(internal_particles[p, m, d+1], x_prev[d+1], S[d+1, d+1]))
                  )
              }
            }
            Z[p] = Z[p] * mean(v)
            
            # On sample des couples selon ces poids
            V <- c(t(v)) # matrice aplatie en vecteur en concaténant les lignes
            sample_indexes <- sample(x=M**2, size = M, replace = T, prob = V)
            
            # %/% (division entière) ; %% (modulo)
            # 12 <-> (2,2) and 20 <-> (2, 10)
            for(i in 1:length(sample_indexes)){
              if(sample_indexes[i]%%M == 0){
                n <- sample_indexes[i]/M
                m <- M
              } else {
                n <- sample_indexes[i]%/%M + 1
                m <- sample_indexes[i]%%M
              }
              updated_internal_particles[p, i, d:(d+1)] <- c(internal_particles[p, n, d], internal_particles[p, m, d+1])
            }
          }
          
          if(k > 1){ # Etape k > 1 : on regroupe par tuple
            internal_particles <- updated_internal_particles
            
            # on regarde toutes les combinaisons possibles de particules et on cacule les poids associés
            v <- array(1, dim=c(M, M))
            deb <- d # debut
            mil <- d+2**(k-1) # milieu
            fin <- d+2**(k)-1 # fin
            for(n in 1:M){
              for(m in 1:M){
                v[n, m] <-  (dmvnorm(c(internal_particles[p, n, deb:(mil-1)], internal_particles[p, m, mil:fin]), x_prev[deb:fin], S[deb:fin, deb:fin]) /
                               (dmvnorm(internal_particles[p, n, deb:(mil-1)], x_prev[deb:(mil-1)], S[deb:(mil-1), deb:(mil-1)]) * 
                               dmvnorm(internal_particles[p, m, mil:fin], x_prev[mil:fin], S[mil:fin, mil:fin])))
              }
            }
            
            Z[p] = Z[p] * mean(v)
            # On sample des couples selon ces poids
            V <- c(t(v)) # matrice en vecteur en concaténant les lignes
            sample_indexes <- sample(x=M**2, size = M, replace = T, prob = V)
            
            # To get couple : index%/%N (division entière) index%%N (modulo)
            
            for(i in 1:length(sample_indexes)){
              if(sample_indexes[i]%%M == 0){
                n <- sample_indexes[i]/M
                m <- M
              }else{
                n <- sample_indexes[i]%/%M + 1
                m <- sample_indexes[i]%%M
              }
              updated_internal_particles[p, i, deb:fin] <- c(internal_particles[p, n, deb:(mil-1)], internal_particles[p, m, mil:fin])
            }
          }
          # Next block
          d <- d+2**k
        }
      }
    }
  }
  
  internal_filter_output = lattice_internal_filter(N, M, D, mu + rho*(particles[, , t-1] - mu) , S, y[, t])
  Z = internal_filter_output$Likelihood
  Z = c(Z/sum(Z))
  updated_internal_particles = internal_filter_output$Particles
  
  resampling_ancestor = sample(N, N, replace=T, Z)
  sampling = sample(M, N, replace=T, matrix(1/M, M))
  for(p in 1:N){
    particles[p, , t] <- updated_internal_particles[sampling[p], , resampling_ancestor[p]]
    #particles[p, , t] <- updated_internal_particles[resampling_ancestor[p], sampling[p], ]
  }
}

for(d in 1:D**2){
  plot(x[d, ])
  lines(x[d, ])
  lines(apply(particles[, d, ], 2, mean), col='red')
}