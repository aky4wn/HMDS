library(MASS)
library(rTensor)
library(abind)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(rstan)
library(reshape2)
library(gridExtra)
library(smacof)
library(car)
library(truncnorm)
library(stringr)



load("Beethoven-hellinger.rda")

N <- dim(hellinger.dist)[1] # Number of orchestras
M <- dim(hellinger.dist)[3] # Number of pieces
K <- dim(hellinger.dist)[4] # Number of metrics
r <- 9 # Dimension of embedding


## Prior values
a <- 0.01
b <- 0.01
alpha <- 1 ## prior for Inv-Gamma for sigma_k^2
beta <- 1
it <- 30000

## Tempo
D <- hellinger.dist[,,,1]/max(hellinger.dist[,,,1]) ## Scale to be between 0-1 
D <- lapply(seq(dim(D)[3]), function(x) D[ , , x])

D.avg <- Reduce("+", D) / length(D)
Y <- cmdscale(D.avg,k = r) 
v <- diag((1/(N-1))*(t(Y) - colMeans(Y))%*%t((t(Y) - colMeans(Y))))

data.tn <- list(N = N, r = r, M = M, D = D, a = a, b = b, alpha = alpha, beta = beta, v = v)
fit.tempo <- stan(file = 'HMDS.stan', data = data.tn, 
                  chains = 1, iter = it, verbose = FALSE, 
                  control = list(max_treedepth = 10,
                                 adapt_delta = 0.93))

## Dynamics
D <- hellinger.dist[,,,2]/max(hellinger.dist[,,,2]) 
D <- lapply(seq(dim(D)[3]), function(x) D[ , , x])

D.avg <- Reduce("+", D) / length(D)
Y <- cmdscale(D.avg,k = r) 
v <- diag((1/(N-1))*(t(Y) - colMeans(Y))%*%t((t(Y) - colMeans(Y))))

data.tn <- list(N = N, r = r, M = M, D = D, a = a, b = b, alpha = alpha, beta = beta, v = v)
fit.volume <- stan(file = 'HMDS.stan', data = data.tn, 
                   chains = 1, iter = it, verbose = FALSE, 
                   control = list(max_treedepth = 10,
                                  adapt_delta = 0.93))


## SF
D <- hellinger.dist[,,,3]/max(hellinger.dist[,,,3]) 
D <- lapply(seq(dim(D)[3]), function(x) D[ , , x])

D.avg <- Reduce("+", D) / length(D)
Y <- cmdscale(D.avg,k = r) 
v <- diag((1/(N-1))*(t(Y) - colMeans(Y))%*%t((t(Y) - colMeans(Y))))

data.tn <- list(N = N, r = r, M = M, D = D, a = a, b = b, alpha = alpha, beta = beta, v = v)
fit.SF <- stan(file = 'HMDS.stan', data = data.tn, 
               chains = 1, iter = it, verbose = FALSE, 
               control = list(max_treedepth = 10,
                              adapt_delta = 0.93))


####################################

# Save stan fit
saveRDS(fit.tempo, "fit-tempo.rds")
saveRDS(fit.volume, "fit-dynamics.rds")
saveRDS(fit.SF, "fit-SF.rds")


####################################