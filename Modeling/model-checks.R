library(MASS)
library(rTensor)
library(abind)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(rstan)
library(reshape2)
library(gridExtra)
library(car)
library(truncnorm)
library(stringr)
library(latex2exp)
set.seed(17)

fit.tempo <- readRDS("fit-tempo.rds")
fit.volume <- readRDS("fit-dynamics.rds")
fit.SF <- readRDS("fit-SF.rds")
load("Beethoven-hellinger.rda")



#########################
### MCMC Diagnostics ####
#########################

## 1. Effective Sample Size
it <- 30000
Neff.tempo <- summary(fit.tempo)$summary[, "n_eff"][1:229] ## only select params

## Only interested in upper triangle of delta_ij matrix
## Values in lower triangle have high ESS since are independent samples (= 0)
inds <- c(2:10, 13:20, 24:30, 35:40, 46:50, 57:60, 68:70, 79:80, 90)
Neff.tempo.delta <- Neff.tempo[130:229][inds]

Neff.tempo.X <- Neff.tempo[1:90]
Neff.tempo.tau <- Neff.tempo[91:127]
Neff.tempo.psi <- Neff.tempo[128]
Neff.tempo.gamma <- Neff.tempo[129]

## Volume
Neff.volume <- summary(fit.volume)$summary[, "n_eff"][1:229] ## only select params

## Only interested in upper triangle of delta_ij matrix
## Values in lower triangle have high ESS since are independent samples (= 0)
inds <- c(2:10, 13:20, 24:30, 35:40, 46:50, 57:60, 68:70, 79:80, 90)
Neff.volume.delta <- Neff.volume[130:229][inds]

Neff.volume.X <- Neff.volume[1:90]
Neff.volume.tau <- Neff.volume[91:127]
Neff.volume.psi <- Neff.volume[128]
Neff.volume.gamma <- Neff.volume[129]

## SF
Neff.SF <- summary(fit.SF)$summary[, "n_eff"][1:229] ## only select params

## Only interested in upper triangle of delta_ij matrix
## Values in lower triangle have high ESS since are independent samples (= 0)
inds <- c(2:10, 13:20, 24:30, 35:40, 46:50, 57:60, 68:70, 79:80, 90)
Neff.SF.delta <- Neff.SF[130:229][inds]

Neff.SF.X <- Neff.SF[1:90]
Neff.SF.tau <- Neff.SF[91:127]
Neff.SF.psi <- Neff.SF[128]
Neff.SF.gamma <- Neff.SF[129]

## Boxplots of Neff

df.delta <- data.frame(Tempo = Neff.tempo.delta, Volume = Neff.volume.delta, 
                       SF = Neff.SF.delta)
df.X <- data.frame(Tempo = Neff.tempo.X, Volume = Neff.volume.X, SF = Neff.SF.X)
df.tau <- data.frame(Tempo = Neff.tempo.tau, Volume = Neff.volume.tau,
                     SF = Neff.SF.tau)


ggplot(melt(df.delta), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(fill = "Metric", x = "Metric", y = "ESS")

ggplot(melt(df.X), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(fill = "Metric", x = "Metric", y = "ESS")

ggplot(melt(df.tau), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  labs(fill = "Metric", x = "Metric", y = "ESS")


round(data.frame(min = apply(df.delta, 2, min),
                 median = apply(df.delta, 2, median),
                 max = apply(df.delta, 2, max)))

round(data.frame(min = apply(df.X, 2, min),
                 median = apply(df.X, 2, median),
                 max = apply(df.X, 2, max)))

round(data.frame(min = apply(df.tau, 2, min),
                 median = apply(df.tau, 2, median),
                 max = apply(df.tau, 2, max)))

round(data.frame(Tempo = Neff.tempo.psi, Volume = Neff.volume.psi, SF = Neff.SF.psi))
round(data.frame(Tempo = Neff.tempo.gamma, Volume = Neff.volume.gamma, SF = Neff.SF.gamma))

## 2. Trace Plots
delta.tempo <- rstan::extract(fit.tempo)$delta
X.post.tempo <- rstan::extract(fit.tempo)$X
tau.tempo <- rstan::extract(fit.tempo)$tau
psi.tempo <- rstan::extract(fit.tempo)$psi
gamma.tempo <- rstan::extract(fit.tempo)$gamma

delta.volume <- rstan::extract(fit.volume)$delta
X.post.volume <- rstan::extract(fit.volume)$X
tau.volume <- rstan::extract(fit.volume)$tau
psi.volume <- rstan::extract(fit.volume)$psi
gamma.volume <- rstan::extract(fit.volume)$gamma

delta.SF <- rstan::extract(fit.SF)$delta
X.post.SF <- rstan::extract(fit.SF)$X
tau.SF <- rstan::extract(fit.SF)$tau
psi.SF <- rstan::extract(fit.SF)$psi
gamma.SF <- rstan::extract(fit.SF)$gamma

## Gamma
gamma.trace <- data.frame(cbind(gamma.tempo, gamma.volume, gamma.SF, 1:15000))
colnames(gamma.trace) <- c("Tempo", "Volume", "SF", "Iteration")
ggplot(melt(gamma.trace, id = c("Iteration")), aes(x = Iteration, y = value, 
                                                   color = variable)) +
  geom_line() +
  theme(text = element_text(size=18)) +
  labs(y = TeX("$\\gamma$"), color = "Metric")


## Psi
psi.trace <- data.frame(cbind(psi.tempo, psi.volume, psi.SF, 1:15000))
colnames(psi.trace) <- c("Tempo", "Volume", "SF", "Iteration")
ggplot(melt(psi.trace, id = c("Iteration")), aes(x = Iteration, y = value, 
                                                 color = variable)) +
  geom_line() +
  theme(text = element_text(size=18)) +
  labs(y = TeX("$\\psi$"), color = "Metric")


## delta - Tempo
it = 30000
delta.post <- matrix(0, nrow = it/2, ncol = 45)
for(i in 1:(it/2)){
  delta.post[i,] <- delta.tempo[i,,][upper.tri(delta.tempo[i,,])]
}

## Order by posterion median
delta.post <- delta.post[,order(apply(delta.post, 2, median))]
## Plot every third entry
delta.plot <- data.frame(delta.post[, seq(1, 45, 5)])

delta.plot$Iteration <- 1:(it/2)

ggplot(melt(delta.plot, id = c("Iteration")), aes(x = Iteration, y = value, 
                                                  color = variable)) +
  geom_line(aes(alpha = 0.8)) +
  labs(y = TeX("$\\delta_{ij}$"), color = "Metric") +
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")


## delta - Volume/Dynamics
delta.post <- matrix(0, nrow = it/2, ncol = 45)
for(i in 1:(it/2)){
  delta.post[i,] <- delta.volume[i,,][upper.tri(delta.volume[i,,])]
}

## Order by posterion median
delta.post <- delta.post[,order(apply(delta.post, 2, median))]
## Plot every third entry
delta.plot <- data.frame(delta.post[, seq(1, 45, 5)])

delta.plot$Iteration <- 1:(it/2)

ggplot(melt(delta.plot, id = c("Iteration")), aes(x = Iteration, y = value, 
                                                  color = variable)) +
  geom_line(aes(alpha = 0.8)) +
  labs(y = TeX("$\\delta_{ij}$"), color = "Metric") +
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")


## delta - SF
delta.post <- matrix(0, nrow = it/2, ncol = 45)
for(i in 1:(it/2)){
  delta.post[i,] <- delta.SF[i,,][upper.tri(delta.SF[i,,])]
}

## Order by posterion median
delta.post <- delta.post[,order(apply(delta.post, 2, median))]
## Plot every third entry
delta.plot <- data.frame(delta.post[, seq(1, 45, 5)])

delta.plot$Iteration <- 1:(it/2)

ggplot(melt(delta.plot, id = c("Iteration")), aes(x = Iteration, y = value, 
                                                  color = variable)) +
  geom_line(aes(alpha = 0.8)) +
  labs(y = TeX("$\\delta_{ij}$"), color = "Metric") +
  theme(text = element_text(size=18)) +
  theme(legend.position = "none")


## Tau
tau.plot <- data.frame(tau.tempo)
colnames(tau.plot) <- dimnames(hellinger.dist)[[3]]
tau.plot$Iteration <- 1:(it/2)
ggplot(melt(tau.plot[,c(1:8, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  theme(text = element_text(size=18)) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(9:16, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(17:25, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(26:38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


## Volume/Dynamics
tau.plot <- data.frame(tau.volume)
colnames(tau.plot) <- dimnames(hellinger.dist)[[3]]
tau.plot$Iteration <- 1:(it/2)
ggplot(melt(tau.plot[,c(1:8, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  theme(text = element_text(size=18)) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(9:16, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(17:25, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(26:38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


## SF
tau.plot <- data.frame(tau.SF)
colnames(tau.plot) <- dimnames(hellinger.dist)[[3]]
tau.plot$Iteration <- 1:(it/2)
ggplot(melt(tau.plot[,c(1:8, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  theme(text = element_text(size=18)) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(9:16, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(17:25, 38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


ggplot(melt(tau.plot[,c(26:38)], id = "Iteration"), 
       aes(x = Iteration, y = value, color = variable)) +
  geom_line(alpha = 0.9) +
  labs(y = TeX("$\\tau_{p}$"), color = "Piece") 


#### 3. Posterior Predictive Checks

## Tempo - Sampling Model

y_star <- array(0, dim = c(15000, 10, 10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,1]/max(hellinger.dist[,,,1])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.tempo, 
                                 rate = psi.tempo/(tau.tempo[,p]*delta.tempo[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

### Posterior pred for one pair across all pieces
post.pred <- data.frame(y_star[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df <-  melt(post.pred)
plot.df$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                      rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                      rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

ggplot(plot.df, aes(x = variable, y = value, fill = Symphony)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), color = 'navy') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\e_{ijp}$"), x = "Piece") +
  theme(text = element_text(size=16)) 

## Posterior pred few pieces all pairs
post.pred <- data.frame(y_star[,1,2:10,34])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")

post.pred <- data.frame(y_star[,1,2:10,1])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")


## Volume/Dynamics - Sampling Model

y_star <- array(0, dim = c(15000, 10, 10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,2]/max(hellinger.dist[,,,2])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.volume, 
                                 rate = psi.volume/(tau.volume[,p]*delta.volume[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

### Posterior pred for one pair across all pieces
post.pred <- data.frame(y_star[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df <-  melt(post.pred)
plot.df$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                      rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                      rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

ggplot(plot.df, aes(x = variable, y = value, fill = Symphony)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), color = 'navy') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\e_{ijp}$"), x = "Piece") +
  theme(text = element_text(size=16)) 

## Posterior pred few pieces all pairs
post.pred <- data.frame(y_star[,1,2:10,34])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")

post.pred <- data.frame(y_star[,1,2:10,1])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")


## SF - Sampling Model

y_star <- array(0, dim = c(15000, 10, 10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,3]/max(hellinger.dist[,,,3])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.SF, 
                                 rate = psi.SF/(tau.SF[,p]*delta.SF[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

### Posterior pred for one pair across all pieces
post.pred <- data.frame(y_star[,1,10,])
colnames(post.pred) <- dimnames(hellinger.dist)[[3]]
plot.df <-  melt(post.pred)
plot.df$Symphony <- c(rep("No1", 15000*4), rep("No2", 15000*4), rep("No3", 15000*4),
                      rep("No4", 15000*4), rep("No5", 15000*4), rep("No6", 15000*5),
                      rep("No7", 15000*4), rep("No8", 15000*4), rep("No9", 15000*4))

ggplot(plot.df, aes(x = variable, y = value, fill = Symphony)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0), color = 'navy') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\e_{ijp}$"), x = "Piece") +
  theme(text = element_text(size=16)) 

## Posterior pred few pieces all pairs
post.pred <- data.frame(y_star[,1,2:10,34])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")

post.pred <- data.frame(y_star[,1,2:10,1])
colnames(post.pred) <- dimnames(hellinger.dist)[[1]][2:10]
ggplot(melt(post.pred), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(y = TeX("$\\y^*_{ijp}$"), x = '') +
  theme(legend.position = "none")

#########################################
####### Evaluate Hierarchical Model ####
########################################

## Tempo

X.del <- array(0, dim = c(15000, 10,10))
for(i in 1:15000){ ## get Euclidean dist based on posteriors of X_i's
  X.del[i,,] <- as.matrix(dist(X.post.tempo[i,,]))
}

## Sample deltas
del <- array(0, dim = c(15000, 10,10))
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      del[,i,j] <- 1/rgamma(15000, shape = gamma.tempo, rate = (gamma.tempo + 1)*X.del[,i,j])
    }
  }
}

## Sample y_ijp's based on the sampled deltas

y_star <- array(0, dim = c(15000, 10,10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,1]/max(hellinger.dist[,,,1])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.tempo, 
                                 rate = psi.tempo/(tau.tempo[,p]*del[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

## Average across pieces
y_star.avg <- apply(y_star, c(1,2,3), mean)

## Take upper triangle
y_star <- matrix(0, nrow = 15000, ncol = 45)
for(i in 1:15000){
  y_star[i,] <- y_star.avg[i,,][upper.tri(y_star.avg[i,,])]
}


ind <- which( upper.tri(y_star.avg[i,,],diag=F) , arr.ind = TRUE )
short.labels <- c("AncientMusic", "Berlin-R", "Berlin-vK", "Chicago",
                  "Leipzig", "LSO", "NBC", "NYPhil", "Philadelphia", "Vienna")
colnames(y_star) <- paste(short.labels[ind[,1]], 
                          short.labels[ind[,2]], sep = "_vs._")


ggplot(melt(y_star), aes(x = Var2, y = value)) +
  geom_boxplot(fill = "#C77CFF") + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = TeX("$\\bar{e}_{ij}$"), x = "Orchestra Pairs") 



## Correlation between ||X_i - X_j|| and obs avg. over pieces
X.del <- apply(X.del, c(2,3), mean)
obs.avg <- apply(D, c(1,2), mean)
plot(X.del[upper.tri(X.del)], obs.avg[upper.tri(obs.avg)])

corr.df <- data.frame(Post.Mean.Tempo = X.del[upper.tri(X.del)],
                      Obs.Mean.Tempo = obs.avg[upper.tri(obs.avg)])


## Volume/Dynamics

X.del <- array(0, dim = c(15000, 10,10))
for(i in 1:15000){ ## get Euclidean dist based on posteriors of X_i's
  X.del[i,,] <- as.matrix(dist(X.post.volume[i,,]))
}

## Sample deltas
del <- array(0, dim = c(15000, 10,10))
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      del[,i,j] <- 1/rgamma(15000, shape = gamma.volume, rate = (gamma.volume + 1)*X.del[,i,j])
    }
  }
}

## Sample y_ijp's based on the sampled deltas

y_star <- array(0, dim = c(15000, 10,10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,2]/max(hellinger.dist[,,,2])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.volume, 
                                 rate = psi.volume/(tau.volume[,p]*del[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

## Average across pieces
y_star.avg <- apply(y_star, c(1,2,3), mean)

## Take upper triangle
y_star <- matrix(0, nrow = 15000, ncol = 45)
for(i in 1:15000){
  y_star[i,] <- y_star.avg[i,,][upper.tri(y_star.avg[i,,])]
}


ind <- which( upper.tri(y_star.avg[i,,],diag=F) , arr.ind = TRUE )
short.labels <- c("AncientMusic", "Berlin-R", "Berlin-vK", "Chicago",
                  "Leipzig", "LSO", "NBC", "NYPhil", "Philadelphia", "Vienna")
colnames(y_star) <- paste(short.labels[ind[,1]], 
                          short.labels[ind[,2]], sep = "_vs._")


ggplot(melt(y_star), aes(x = Var2, y = value)) +
  geom_boxplot(fill = "#00BFC4") + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = TeX("$\\bar{e}_{ij}$"), x = "Orchestra Pairs")

## Correlation between ||X_i - X_j|| and obs avg. over pieces
X.del <- apply(X.del, c(2,3), mean)
obs.avg <- apply(D, c(1,2), mean)
plot(X.del[upper.tri(X.del)], obs.avg[upper.tri(obs.avg)])

corr.df$Post.Mean.Volume <- X.del[upper.tri(X.del)]
corr.df$Obs.Mean.Volume <- obs.avg[upper.tri(obs.avg)]


## SF

X.del <- array(0, dim = c(15000, 10,10))
for(i in 1:15000){ ## get Euclidean dist based on posteriors of X_i's
  X.del[i,,] <- as.matrix(dist(X.post.SF[i,,]))
}

## Sample deltas
del <- array(0, dim = c(15000, 10,10))
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      del[,i,j] <- 1/rgamma(15000, shape = gamma.SF, rate = (gamma.SF + 1)*X.del[,i,j])
    }
  }
}

## Sample y_ijp's based on the sampled deltas

y_star <- array(0, dim = c(15000, 10,10, 37))
## Scale true data between 0 and 1
D <- hellinger.dist[,,,3]/max(hellinger.dist[,,,3])
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        ## Subtract off obs data
        y_star[,i,j,p] <- rgamma(15000, shape = psi.SF, 
                                 rate = psi.SF/(tau.SF[,p]*del[,i,j])) - D[i,j,p]
      }
    }
  }
}

q25 <- apply(y_star, c(2,3,4), quantile, 0.025)
q975 <- apply(y_star, c(2,3,4), quantile, 0.975)
## Check if 0 in quantile
count <- 0
for(i in 1:10){
  for(j in 1:10){
    if(j > i){
      for(p in 1:37){
        if(q25[i,j,p] <0 & q975[i,j,p] > 0){
          count <- count + 1
        }
      }
    }
  }
}

## Average across pieces
y_star.avg <- apply(y_star, c(1,2,3), mean)

## Take upper triangle
y_star <- matrix(0, nrow = 15000, ncol = 45)
for(i in 1:15000){
  y_star[i,] <- y_star.avg[i,,][upper.tri(y_star.avg[i,,])]
}


ind <- which( upper.tri(y_star.avg[i,,],diag=F) , arr.ind = TRUE )
short.labels <- c("AncientMusic", "Berlin-R", "Berlin-vK", "Chicago",
                  "Leipzig", "LSO", "NBC", "NYPhil", "Philadelphia", "Vienna")
colnames(y_star) <- paste(short.labels[ind[,1]], 
                          short.labels[ind[,2]], sep = "_vs._")


ggplot(melt(y_star), aes(x = Var2, y = value)) +
  geom_boxplot(fill = "#F8765D") + 
  geom_hline(aes(yintercept = 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(y = TeX("$\\bar{e}_{ij}$"), x = "Orchestra Pairs")

 
## Correlation between ||X_i - X_j|| and obs avg. over pieces
X.del <- apply(X.del, c(2,3), mean)
obs.avg <- apply(D, c(1,2), mean)
plot(X.del[upper.tri(X.del)], obs.avg[upper.tri(obs.avg)])

corr.df$Post.Mean.SF <- X.del[upper.tri(X.del)]
corr.df$Obs.Mean.SF <- obs.avg[upper.tri(obs.avg)]

corr.df.obs <- corr.df[,c(2,4,6)]
colnames(corr.df.obs) <- c("Tempo", "Volume", "SF")
corr.df.post <- corr.df[,c(1,3,5)]
colnames(corr.df.post) <- c("Tempo", "Volume", "SF")

plot.obs <- melt(corr.df.obs)
colnames(plot.obs) <- c("Metric", "Observed")

plot.post <- melt(corr.df.post)
colnames(plot.post) <- c("Metric", "Posterior")

plot.df <- cbind(plot.post, plot.obs)[,-3]

ggplot(plot.df, aes(x = Observed, y = Posterior, color = Metric, shape = Metric)) +
  geom_point() +
  labs(x = TeX("Observed: $\\bar{D}_{ij}$"),
       y = TeX("Posterior: $\\bar{||X_i - X_j||_2}$"))
