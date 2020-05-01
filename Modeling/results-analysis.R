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
library(RColorBrewer)
library(latex2exp)
set.seed(17)

fit.tempo <- readRDS("fit-tempo.rds")
fit.volume <- readRDS("fit-dynamics.rds")
fit.SF <- readRDS("fit-SF.rds")
load("Beethoven-hellinger.rda")


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


####################################

## Look at psi
psi.df <- data.frame(tempo = c(psi.tempo), dynamics = c(psi.volume), SF = c(psi.SF),
                     it = 1:it)

ggplot(melt(psi.df, id = "it"), aes(x = variable, y = log(value), fill = variable)) +
  geom_violin() +
  labs(x = "Metric", y = TeX("$\\log(\\psi)$"), fill = "Metric")

## Gamma
gamma.df <- data.frame(tempo = c(gamma.tempo), dynamics = c(gamma.volume), SF = c(gamma.SF),
                       it = 1:it)

ggplot(melt(gamma.df, id = "it"), aes(x = variable, y = log(value), fill = variable)) +
  geom_violin() +
  labs(x = "Metric", y = TeX("$\\log(\\gamma)$"), fill = "Metric")

## Tau

tau1.post <- cbind(tau.tempo[,34], tau.volume[,34], tau.SF[,34])
tau2.post <- cbind(tau.tempo[,35], tau.volume[,35], tau.SF[,35])
tau3.post <- cbind(tau.tempo[,36], tau.volume[,36], tau.SF[,36])
tau4.post <- cbind(tau.tempo[,37], tau.volume[,37], tau.SF[,37])


colnames(tau1.post) <- c("Tempo", "Dynamics", "SF")
colnames(tau2.post) <- c("Tempo", "Dynamics", "SF")
colnames(tau3.post) <- c("Tempo", "Dynamics", "SF")
colnames(tau4.post) <- c("Tempo", "Dynamics", "SF")


p1 <- ggplot(melt(tau1.post), aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  labs(x = "Metric", y = TeX("$\\tau_p$"),  fill = "Metric") +
  ggtitle("No9-01")

p2 <- ggplot(melt(tau2.post), aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  labs(x = "Metric", y = TeX("$\\tau_p$"),  fill = "Metric")+
  ggtitle("No9-02")

p3 <- ggplot(melt(tau3.post), aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  labs(x = "Metric", y = TeX("$\\tau_p$"),  fill = "Metric")+
  ggtitle("No9-03")

p4 <- ggplot(melt(tau4.post), aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  labs(x = "Metric", y = TeX("$\\tau_p$"),  fill = "Metric")+
  ggtitle("No9-04")


p <- grid.arrange(p1, p2, p3, p4,   nrow = 2)



## Tau plot by metric
tau.df <- data.frame(tau.volume)
colnames(tau.df) <- dimnames(hellinger.dist)[[3]]

tau.df <- melt(tau.df)
sp.df <- str_split_fixed(tau.df$variable, "-", 2)
colnames(sp.df) <- c("Symphony", "Movement")
tau.df <- cbind(tau.df, sp.df)

ggplot(tau.df, aes(x = Movement, y = value, fill = Movement)) +
  geom_boxplot() +
  facet_wrap(~Symphony) +
  labs(y = TeX("$\\tau$")) +
  theme(text = element_text(size=18))



############################################
## Plot delta posterior values

value.true <- apply(delta.tempo, c(2,3), mean) ## posterior mean

plot.df <- delta.tempo
dimnames(plot.df)[[2]] <- dimnames(hellinger.dist)[[1]]
dimnames(plot.df)[[3]] <- dimnames(hellinger.dist)[[1]]

plots <- list()
count <- 1
col.order <- order(value.true[upper.tri(value.true)])
colfunc <- colorRampPalette(c("#132B43", "#56B1F7"))
cols <- colfunc(45)
for(j in 2:N){
  for(i in 1:(j-1)){
    ind <- which(col.order == count)
    
    ## Titles
    if(count == 1){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = TeX("$\\delta_{ij}$")) +
        #ggtitle("Berlin-Rattle") +
        theme(plot.title = element_text(size=20)) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
              text=element_text(size=20))
    }
    else if(count == 2){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("Berlin-vK") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
              text=element_text(size=20))
    }
    else if(count == 4){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("Chicago") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 7){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("Leipzig") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 11){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("LSO") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 16){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("NBC") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 22){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("NYPhil") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 29){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("Philadelphia") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count == 37){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "", y = "") +
        #ggtitle("Vienna") +
        theme(plot.title = element_text(size=20))+
        theme(axis.text.x=element_text(angle=45, hjust=1),
              text=element_text(size=20))
    }
    else if(count %in% c(3,6,10,15,21,28,36,45)){
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = TeX("$\\delta_{ij}$")) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
              text=element_text(size=20))
    }
    else{
      p <- ggplot(data.frame(value = plot.df[,i,j]), aes(x = value)) +
        geom_histogram(fill = cols[ind]) +
        geom_vline(xintercept = value.true[i,j]) +
        #ylim(0, 2500) +
        xlim(0,0.05) +
        labs(x = "")+
        theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
              text=element_text(size=20))
    }
    plots[[count]] <- p
    count <- count + 1
  }
}


m <- matrix(NA, 10, 10)
m[upper.tri(m, diag = F)] <- 1:45
mat <- grid.arrange(grobs = plots, layout_matrix = m)


## Point estimates - using posterior means for delta
mean.tempo <- apply(delta.tempo, c(2,3), mean)
mean.volume <- apply(delta.volume, c(2,3), mean)
mean.SF <- apply(delta.SF, c(2,3), mean)

## Hierarchical Clustering on delta
# Dissimilarity matrix
delta_avg <- mean.volume
delta_avg[lower.tri(delta_avg)]  <- t(delta_avg)[lower.tri(delta_avg)]
# Hierarchical clustering using Complete Linkage
colnames(delta_avg) <- dimnames(hellinger.dist)[[1]]
hc1 <- hclust(as.dist(delta_avg), method = "complete" )




## Heatplots of delta_avg

df <- data.frame(mean.volume)
df[lower.tri(df)] <-NA
diag(df) <- NA
colnames(df) <- dimnames(hellinger.dist)[[1]]
df$Orchestra <- factor(dimnames(hellinger.dist)[[1]], 
                       levels = rev(dimnames(hellinger.dist)[[1]]))

ggplot(melt(df), aes(x = variable, y = Orchestra, fill= value)) + 
  geom_tile(colour = "white") +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE),
        text = element_text(size=18)) +
  labs(x = "", y = "", fill = TeX("$\\bar{\\delta}_{ij}$"))+
  scale_fill_distiller(palette = "RdPu", direction=-1)




## Plot ||X_i - X_j||

X.del <- array(0, dim = c(15000, 10,10))
for(i in 1:15000){ ## get Euclidean dist based on posteriors of X_i's
  X.del[i,,] <- as.matrix(dist(X.post.SF[i,,]))
}

mean.del <- apply(X.del, c(2,3), mean)
df <- data.frame(mean.del)
df[lower.tri(df)] <-NA
diag(df) <- NA
colnames(df) <- dimnames(hellinger.dist)[[1]]
df$Orchestra <- factor(dimnames(hellinger.dist)[[1]], 
                       levels = rev(dimnames(hellinger.dist)[[1]]))

ggplot(melt(df), aes(x = variable, y = Orchestra, fill= value)) + 
  geom_tile(colour = "white") +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE),
        text = element_text(size=18)) +
  labs(x = "", y = "", fill = TeX("$||\\bar{X_i - X_j}||$"))+
  scale_fill_distiller(palette = "RdPu", direction=-1)




## Boxplots of similar/dissimilar pairs
tempo.df <- data.frame("Vienna_BerlinRattle" = delta.tempo[,2,10],
                       "Chicago_Philadelphia" = delta.tempo[,4,9],
                       "BerlinRattle_BerlinvonKarajan" = delta.tempo[,2,3],
                       "BerlinvonKarajan_NY" = delta.tempo[,3,8])
plot.df <- melt(tempo.df)
plot.df$Similarity <- rep(c("Similar", "Different"),
                          each = 30000)
ggplot(plot.df, aes(x = variable, y = value, fill = Similarity)) +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE)) +
  geom_boxplot()+
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$"))

## Volume
df <- data.frame("Chicago_NY" = delta.volume[,4,8],
                 "AncientMusic_LSO" = delta.volume[,1,6])
plot.df <- melt(df)
plot.df$Similarity <- rep(c("Similar", "Different"),
                          each = 15000)
p1 <- ggplot(plot.df, aes(x = variable, y = value, fill = Similarity)) +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE)) +
  geom_boxplot()+
  #facet_wrap(~Metric) +
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$")) +
  ggtitle("Volume")

## SF
df <- data.frame("BerlinvonKarajan_Leipzig" = delta.SF[,3,5],
                 "NBC_Vienna" = delta.SF[,7,10])
plot.df <- melt(df)
plot.df$Similarity <- rep(c("Similar", "Different"),
                          each = 15000)
p2 <- ggplot(plot.df, aes(x = variable, y = value, fill = Similarity)) +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE)) +
  geom_boxplot()+
  #facet_wrap(~Metric) +
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$")) +
  ggtitle("Spectral Flatness")

## spec5
df <- data.frame("BerlinRattle_Vienna" = delta.spec5[,2,10],
                 "LSO_NBC" = delta.spec5[,6,7])
plot.df <- melt(df)
plot.df$Similarity <- rep(c("Similar", "Different"),
                          each = 15000)
p3 <- ggplot(plot.df, aes(x = variable, y = value, fill = Similarity)) +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE)) +
  geom_boxplot()+
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$")) +
  ggtitle("Volume - A4-A5")
## SF5
df <- data.frame("BerlinRattle_Vienna" = delta.SF5[,2,10],
                 "AncientMusic_NBC" = delta.SF5[,1,7])
plot.df <- melt(df)
plot.df$Similarity <- rep(c("Similar", "Different"),
                          each = 15000)
p4 <- ggplot(plot.df, aes(x = variable, y = value, fill = Similarity)) +
  theme(axis.text.x = element_text(angle=45, hjust=TRUE)) +
  geom_boxplot()+
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$")) +
  ggtitle("Spectral Flatness - A4-A5")
p <- grid.arrange(p1, p2, p3, p4, nrow = 2)



#################################
## Embed X_i - X_j 
################################


## Tempo

library(Rtsne)

## Align X_i to overall average
D <- hellinger.dist[,,,1]/max(hellinger.dist[,,,1])
X.target <- cmdscale(apply(D, c(1,2), mean), k = 9)
## Align via Procrustes to mean
X.embed <- array(0, dim = c(15000, 10,9))
for(i in 1:15000){
  X.embed[i,,] <- Procrustes(X.target, X.post.tempo[i,,])$Yhat
}


N <- 10
X.plot <- data.frame(array(X.post.tempo[c(TRUE, FALSE, FALSE),,], dim = c(15000/3*N, 9)))
X.plot$Orchestra <- rep(dimnames(hellinger.dist)[[1]], each = 15000/3)

tsne <- Rtsne(X.plot[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
df.plot <- data.frame(tsne$Y)
df.plot$Orchestra <- X.plot$Orchestra

X.mean <- aggregate(df.plot[, 1:2], list(df.plot$Orchestra), mean)

ggplot(melt(df.plot, id = c("X1", "X2"))) +
  geom_point(aes(x = X1, y = X2, group = value, color = value, shape = value),
             size = 2.5) +
  scale_shape_manual(values=c("0", "1", "2", "3", "4", "5", "6", "7", "8","9"),
                     name = "Orchestra") +
  geom_text(data = data.frame(X.mean), aes(x = X1, y = X2),
            colour="black", size=4, 
            label = c("0", "1", "2", "3", "4", "5", "6", "7", "8","9")) +
  labs(color = "Orchestra", x = "X1", y = "X2") 

## Data really needs to be represented in > 2 dim space (X vectors), tsne helps visualize
## since nonlinear


## Volume

## Align X_i to overall average
D <- hellinger.dist[,,,2]/max(hellinger.dist[,,,2])
X.target <- cmdscale(apply(D, c(1,2), mean), k = 9)
## Align via Procrustes to mean
X.embed <- array(0, dim = c(15000, 10,9))
for(i in 1:15000){
  X.embed[i,,] <- Procrustes(X.target, X.post.volume[i,,])$Yhat
}


N <- 10
X.plot <- data.frame(array(X.post.tempo[c(TRUE, FALSE, FALSE),,], dim = c(15000/3*N, 9)))
X.plot$Orchestra <- rep(dimnames(hellinger.dist)[[1]], each = 15000/3)

tsne <- Rtsne(X.plot[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
df.plot <- data.frame(tsne$Y)
df.plot$Orchestra <- X.plot$Orchestra

X.mean <- aggregate(df.plot[, 1:2], list(df.plot$Orchestra), mean)

ggplot(melt(df.plot, id = c("X1", "X2"))) +
  geom_point(aes(x = X1, y = X2, group = value, color = value, shape = value),
             size = 2.5) +
  scale_shape_manual(values=c("0", "1", "2", "3", "4", "5", "6", "7", "8","9"),
                     name = "Orchestra") +
  geom_text(data = data.frame(X.mean), aes(x = X1, y = X2),
            colour="black", size=4, 
            label = c("0", "1", "2", "3", "4", "5", "6", "7", "8","9")) +
  labs(color = "Orchestra", x = "X1", y = "X2") 

## SF

## Align X_i to overall average
D <- hellinger.dist[,,,3]/max(hellinger.dist[,,,3])
X.target <- cmdscale(apply(D, c(1,2), mean), k = 9)
## Align via Procrustes to mean
X.embed <- array(0, dim = c(15000, 10,9))
for(i in 1:15000){
  X.embed[i,,] <- Procrustes(X.target, X.post.SF[i,,])$Yhat
}


N <- 10
X.plot <- data.frame(array(X.post.tempo[c(TRUE, FALSE, FALSE),,], dim = c(15000/3*N, 9)))
X.plot$Orchestra <- rep(dimnames(hellinger.dist)[[1]], each = 15000/3)

tsne <- Rtsne(X.plot[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
df.plot <- data.frame(tsne$Y)
df.plot$Orchestra <- X.plot$Orchestra

X.mean <- aggregate(df.plot[, 1:2], list(df.plot$Orchestra), mean)

ggplot(melt(df.plot, id = c("X1", "X2"))) +
  geom_point(aes(x = X1, y = X2, group = value, color = value, shape = value),
             size = 2.5) +
  scale_shape_manual(values=c("0", "1", "2", "3", "4", "5", "6", "7", "8","9"),
                     name = "Orchestra") +
  geom_text(data = data.frame(X.mean), aes(x = X1, y = X2),
            colour="black", size=4, 
            label = c("0", "1", "2", "3", "4", "5", "6", "7", "8","9")) +
  labs(color = "Orchestra", x = "X1", y = "X2") 

## Boxplots for delta_ij
vienna.tempo <- data.frame(delta.tempo[,1:9, 10])
vienna.volume <- data.frame(delta.volume[,1:9, 10])
vienna.SF <- data.frame(delta.SF[,1:9, 10])

colnames(vienna.tempo) <- dimnames(hellinger.dist)[[1]][1:9]
colnames(vienna.volume) <- dimnames(hellinger.dist)[[1]][1:9]
colnames(vienna.SF) <- dimnames(hellinger.dist)[[1]][1:9]

vienna.tempo$Metric <- 'Tempo'
vienna.volume$Metric <- 'Dynamics'
vienna.SF$Metric <- 'Timbre'

vienna.df <- rbind(vienna.tempo, vienna.volume, vienna.SF)
vienna.df$Metric <- factor(vienna.df$Metric, 
                           levels = c("Tempo", "Dynamics", "Timbre"))
ggplot(melt(vienna.df, id = "Metric"), aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  facet_wrap(~Metric, nrow = 1) +
  labs(x = "Orchestra Pair", y = TeX("$\\delta_{ij}$"), fill = "Orchestra") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        text = element_text(size=18)) 

#### Check Triangle Inequality
tempo.tri <- rep(0, 239)
count <- 0
for(i in 1:9){
  for(j in (i+1):10){
    for(k in (i+1):10){
      if(k != j){
        print(c(i,j, k)) 
        if(k < j){
          tempo.tri[count] <- length(which(delta.tempo[,i,j] <= delta.tempo[,i,k] + 
                                             delta.tempo[, k,j]))
        }
        else{
          tempo.tri[count] <- length(which(delta.tempo[,i,j] <= delta.tempo[,i,k] + 
                                             delta.tempo[, j,k])) 
        }
        count <- count + 1
      }
    }
  }
}

all(tempo.tri == 15000)


vol.tri <- rep(0, 239)
count <- 0
for(i in 1:9){
  for(j in (i+1):10){
    for(k in (i+1):10){
      if(k != j){
        print(c(i,j, k)) 
        if(k < j){
          vol.tri[count] <- length(which(delta.volume[,i,j] <= delta.volume[,i,k] + 
                                           delta.volume[, k,j]))
        }
        else{
          vol.tri[count] <- length(which(delta.volume[,i,j] <= delta.volume[,i,k] + 
                                           delta.volume[, j,k])) 
        }
        count <- count + 1
      }
    }
  }
}

all(vol.tri == 15000)


sf.tri <- rep(0, 239)
count <- 0
for(i in 1:9){
  for(j in (i+1):10){
    for(k in (i+1):10){
      if(k != j){
        print(c(i,j, k)) 
        if(k < j){
          sf.tri[count] <- length(which(delta.SF[,i,j] <= delta.SF[,i,k] + 
                                          delta.SF[, k,j]))
        }
        else{
          sf.tri[count] <- length(which(delta.SF[,i,j] <= delta.SF[,i,k] + 
                                          delta.SF[, j,k])) 
        }
        count <- count + 1
      }
    }
  }
}

all(sf.tri == 15000)
