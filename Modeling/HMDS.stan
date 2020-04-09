//  BMDS - w/reps and Separate Euclidean Assumption

data {
  int<lower=0> N; // Number of orchestras
  int<lower=0> r; // Dimension of embedding
  int<lower=0> M; // number of pieces
  matrix<lower=0>[N, N] D[M]; // observed distances, list by piece
  real<lower=0> a; // prior for Gamma for psi and gamma
  real<lower=0> b; // prior for Gamma for psi and gamma
  real<lower=0> alpha; // hyper-prior for Inv-Gamma for tau
  real<lower=0> beta; // hyper-prior for Inv-Gamma for tau
  vector<lower=0>[r] v; // covariance matrix for X, vI_r
}

parameters {
  matrix[N, r] X; // embedding vectors
  vector<lower=0>[M] tau; // scaling for amount of variation for each piece
  real<lower=0> psi; // noise variance
  real<lower=0> gamma; // Euclidean assumption
  matrix<lower=0>[N, N] delta; 
  
}
transformed parameters {
  matrix<lower=0>[N, N] dx; //  distance matrix from embedding, orchestras x orchestras

  // Calculate Euclidean Distance
  for(i in 1:N){
    for(j in 1:N){
      if(j > i){
        dx[i, j] = distance(X[i,], X[j,]); //Euclidean distance
      } else{
        dx[i,j] = 0;
      }
    }
  }

}

model {
  
  // likelihood
    for(p in 1:M){
      for(i in 1:N){
        for(j in 1:N){
          if(j>i){
            D[p][i,j] ~ gamma(psi, psi/(tau[p]*delta[i,j]));
            
          } 
        }
      }
    }
    
    // Delta values
    for(i in 1:N){
      for(j in 1:N){
        if(j>i){
          delta[i,j] ~ inv_gamma(gamma, (gamma+1)*dx[i,j]);
        } else{
          delta[i,j] ~ normal(0, 0.000000001); // 0 otherwise
        }
      }
    }
    
   
  // priors
  for(i in 1:N){
    for(d in 1:r){
    X[i, d] ~ normal(0, v[d]); //Lambda is a diagonal matrix
    }
  }


  psi ~ gamma(a, b); 
  gamma ~ gamma(a, b);
  
  
 for(p in 1:M){
    tau[p] ~ inv_gamma(alpha, beta);
  }
  
}

