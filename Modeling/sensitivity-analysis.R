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
library(latex2exp)



load("Beethoven-hellinger.rda")

N <- dim(hellinger.dist)[1] # Number of orchestras
M <- dim(hellinger.dist)[3] # Number of pieces
K <- dim(hellinger.dist)[4] # Number of metrics
r <- 9 # Dimension of embedding

## (1) ######################


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


tau.tempo1 <- rstan::extract(fit.tempo)$tau
psi.tempo1 <- rstan::extract(fit.tempo)$psi
gamma.tempo1 <- rstan::extract(fit.tempo)$gamma
delta.tempo1 <- rstan::extract(fit.tempo)$delta

saveRDS(fit.tempo, "sa-tempo-1.rds")

## (2) ######################


## Prior values
a <- 0.01
b <- 0.01
alpha <- 0.01 ## prior for Inv-Gamma for sigma_k^2
beta <- 0.01
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


tau.tempo2 <- rstan::extract(fit.tempo)$tau
psi.tempo2 <- rstan::extract(fit.tempo)$psi
gamma.tempo2 <- rstan::extract(fit.tempo)$gamma
delta.tempo2 <- rstan::extract(fit.tempo)$delta
saveRDS(fit.tempo, "sa-tempo-2.rds")
## (3) ######################
## Centered about 1, much shorter tails

## Prior values
a <- 0.01
b <- 0.01
alpha <- 10 ## prior for Inv-Gamma for sigma_k^2
beta <- 9
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


tau.tempo3 <- rstan::extract(fit.tempo)$tau
psi.tempo3 <- rstan::extract(fit.tempo)$psi
gamma.tempo3 <- rstan::extract(fit.tempo)$gamma
delta.tempo3 <- rstan::extract(fit.tempo)$delta
saveRDS(fit.tempo, "sa-tempo-3.rds")

## (4) ######################


## Prior values
a <- 1
b <- 1
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


tau.tempo4 <- rstan::extract(fit.tempo)$tau
psi.tempo4 <- rstan::extract(fit.tempo)$psi
gamma.tempo4 <- rstan::extract(fit.tempo)$gamma
delta.tempo4 <- rstan::extract(fit.tempo)$delta

saveRDS(fit.tempo, "sa-tempo-4.rds")


################################################################
### Posterior Predictive Checks to Make Sure Values Constant ###
################################################################

## Eval sampling model

y_star1 <- array(0, dim = c(15000, 10, 10, 37))
y_star2 <- array(0, dim = c(15000, 10, 10, 37))
y_star3 <- array(0, dim = c(15000, 10, 10, 37))
y_star4 <- array(0, dim = c(15000, 10, 10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,1]/max(hellinger.dist[,,,1])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Divide by obs data
        y_star1[,i,j,p] <- rgamma(15000, shape = psi.tempo1, 
                                 rate = psi.tempo1/(tau.tempo1[,p]*delta.tempo1[,i,j]))/D[i,j,p]
        y_star2[,i,j,p] <- rgamma(15000, shape = psi.tempo2, 
                                  rate = psi.tempo2/(tau.tempo2[,p]*delta.tempo2[,i,j]))/D[i,j,p]
        y_star3[,i,j,p] <- rgamma(15000, shape = psi.tempo3, 
                                  rate = psi.tempo3/(tau.tempo3[,p]*delta.tempo3[,i,j]))/D[i,j,p]
        y_star4[,i,j,p] <- rgamma(15000, shape = psi.tempo4, 
                                  rate = psi.tempo4/(tau.tempo4[,p]*delta.tempo4[,i,j]))/D[i,j,p]
      }
    }
  }
}

q25_1 <- apply(y_star1, c(2,3,4), quantile, 0.025)
q975_1 <- apply(y_star1, c(2,3,4), quantile, 0.975)
q25_2 <- apply(y_star2, c(2,3,4), quantile, 0.025)
q975_2 <- apply(y_star2, c(2,3,4), quantile, 0.975)
q25_3 <- apply(y_star3, c(2,3,4), quantile, 0.025)
q975_3 <- apply(y_star3, c(2,3,4), quantile, 0.975)
q25_4 <- apply(y_star4, c(2,3,4), quantile, 0.025)
q975_4 <- apply(y_star4, c(2,3,4), quantile, 0.975)
## Check if 1 in quantile
count1 <- 0
count2 <- 0
count3 <- 0
count4 <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25_1[i,j,p] <1 & q975_1[i,j,p] > 1){
          count1 <- count1 + 1
        }
        if(q25_2[i,j,p] <1 & q975_2[i,j,p] > 1){
          count2 <- count2 + 1
        }
        if(q25_3[i,j,p] <1 & q975_3[i,j,p] > 1){
          count3 <- count3 + 1
        }
        if(q25_4[i,j,p] <1 & q975_4[i,j,p] > 1){
          count4 <- count4 + 1
        }
      }
    }
  }
}

round(c(count1/(37*45), count2/(37*45), count3/(37*45), count4/(37*45)), 2)








post.pred <- data.frame(y_star1[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df1 <-  melt(post.pred)
plot.df1$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                      rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                      rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

plot.df1$Model <- '1'

post.pred <- data.frame(y_star2[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df2 <-  melt(post.pred)
plot.df2$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                       rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                       rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

plot.df2$Model <- '2'

post.pred <- data.frame(y_star3[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df3 <-  melt(post.pred)
plot.df3$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                       rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                       rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

plot.df3$Model <- '3'

post.pred <- data.frame(y_star4[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df4 <-  melt(post.pred)
plot.df4$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                       rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                       rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

plot.df4$Model <- '4'

plot.df <- rbind(plot.df1, plot.df2, plot.df3, plot.df4)
sp.df <- str_split_fixed(plot.df$variable, "-", 2)
colnames(sp.df) <- c("Symphony2", "Movement")
plot.df <- cbind(plot.df, sp.df)


ggplot(plot.df, aes(x = Movement, y = log(value), fill = Model)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), color = 'navy') +
  facet_wrap(~Symphony) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\r_{ijp}$"), x = "Piece") +
  theme(text = element_text(size=16, family="Times")) 
ggsave(filename = "sensitivity.png", height = 20, width = 25, units = "cm")


##### Table of tau_p values

t1.mean <- colMeans(tau.tempo1/max(tau.tempo1))
t1.sd <- apply(tau.tempo1/max(tau.tempo1), 2, sd)
t2.mean <- colMeans(tau.tempo2/max(tau.tempo2))
t2.sd <- apply(tau.tempo2/max(tau.tempo2), 2, sd)
t3.mean <- colMeans(tau.tempo3/max(tau.tempo3))
t3.sd <- apply(tau.tempo3/max(tau.tempo3), 2, sd)
t4.mean <- colMeans(tau.tempo4/max(tau.tempo4))
t4.sd <- apply(tau.tempo4/max(tau.tempo4), 2, sd)

round(t1.mean[1:8],2)
round(t1.sd[1:8],2)
round(t2.mean[1:8],2)
round(t2.sd[1:8],2)
round(t3.mean[1:8],2)
round(t3.sd[1:8],2)
round(t4.mean[1:8],2)
round(t4.sd[1:8],2)


