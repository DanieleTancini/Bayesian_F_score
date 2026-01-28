#In this script we analyze the real data application
#In the first block is considered the Bayesian HMM model, including the multimodel framework and the model averaging. The frequentist logistic model is considered in the same block
#In the second block, we consider the Bayesian logistic model
#In the last block we report the results

#----------------------
#Bayesian HMM and frequentist logistic model
#----------------------

#Set your WD
#setwd("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score")

Y = readRDS("Y_matrix_missing")
X = readRDS("X_complete_new")

#library

library("LaplacesDemon")
library("abind")

#Functions
log_prior_pi = function(x, a, b){
  dbeta(x, a, b, log = TRUE)
}

log_jac = function(phin, phi){
  log(exp(phin)/(1+exp(phin))^2) - log(exp(phi)/(1+exp(phi))^2)
}

log_beta = function(beta, mu, sigma, sigma2, pi_beta_){
  log(pi_beta_*dlaplace(beta, mu, sigma) + (1-pi_beta_)*dlaplace(beta, mu, sigma2))
}

log_alpha = function(alpha, mu, sigma){
  sum(dnorm(alpha, mu, sigma, log = TRUE))
}

log_omega = function(omega, shape, scale){
  sum(dgamma(omega, shape = shape, scale = scale, log = TRUE))
}

log_omega2 = function(ome, vec){
  ome = as.matrix(ome)
  val = 0
  for(i in 1:nrow(ome)){
    for(j in 1:nrow(ome)){
      val = val + dgamma(ome[i,j], shape = vec[i,j], 1, log = TRUE)
    }
  }
  val
}

log_lik = function(y, u, x, alpha, beta, epsilon){
  N = as.numeric(dim(x)[1])
  Ti = as.numeric(dim(x)[3])
  l = 0
  for(i in 1:N){
    for(t in 1:Ti){
      psi = alpha[u[i,t]] + beta%*%x[i,,t]
      p = exp(psi)/(1+exp(psi))
      if(is.nan(p)==TRUE){
        p = 1
      }
      if(p == 0){
        p = p + epsilon
      }
      if(p == 1){
        p = p - epsilon
      }
      l = l + dbinom(y[i,t], 1, p, log = TRUE)
    }
  }
  l
}

log_latent = function(U_, pi, Pi, K){
  ind1 = as.vector(table(factor(U_[,1], levels = 1:K)))
  part1 = sum(ind1*log(pi))
  Ti = ncol(U)
  N = nrow(U)
  ind2 = 0
  for(i in 1:N){
    ind2 = ind2 + as.vector(table(factor(U_[i,1:Ti-1], levels = 1:K), factor(U_[i,2:Ti], levels = 1:K)))
  }
  part2 = sum(ind2*log(as.vector(Pi)))
  part1 + part2
}

Jacobian_omega = function(x){
  sum(log(x))
}

Jacobian_Omega = function(x){
  sum(as.vector(log(x)))
}

Jacobian_split = function(oms, alps, omhs, omsh, omss, dh, e, f, g){
  log(abs(oms*4*alps*prod(omhs)*prod(2*omsh/dh)*4*(omss^3)*e*(e-1)/(f*g)))
}

Jacobian_split_1 = function(oms, alps, omss, e, f, g){
  log(abs(oms*4*alps*4*(omss^3)*e*(e-1)/(f*g)))
}


#Number of classes
K_min = 1
K_max = 10

#Number of iterations
R = 50000

#Dataset defined as N,P,T
#i = 1,...,N and t = 1,...,T

#Training and test set
T_fin = 35
T_m1 = 5
T_f = T_m1+1

for(ii in 1:29){
  Y = readRDS("Y_matrix_missing")
  X = readRDS("X_complete_new")
  
  #same scale
  X[,7,] = 100*X[,7,]
  
  #Windows
  Y_test = Y[,(T_f+ii):T_fin]
  X_test = X[,,(T_f+ii):T_fin]
  Y = Y[,ii:(T_m1+ii)]
  X = X[,,ii:(T_m1+ii)]
  N = as.numeric(nrow(Y))
  Ti = as.numeric(ncol(Y))
  
  K = 10
  
  Y_impute = matrix(0,N,Ti)
  for(t in 1:(T_m1+1)){
    Y_impute[which(is.na(Y[,t])),t] = 1  
  }
  
  #Number of predictors
  P = as.numeric(dim(X)[2])
  
  #K-state
  K_state = rep(NA, R)
  K_state[1] = K
  
  #seed
  set.seed(123456)
  
  #Storage
  Alpha = matrix(NA, K_max, R)
  Beta = matrix(NA, P, R)
  omega_ = matrix(NA, K_max, R)
  Omega_ = array(NA, dim = c(K_max, K_max, R))
  Pi_ = array(NA, dim = c(K_max, K_max, R))
  pi_ = matrix(NA, K_max, R)
  U = array(NA, dim = c(N,Ti,R))
  Pi_beta = matrix(NA, P, R)
  
  #Initial values
  Alpha[1:K,1] = alpha = alphan = rnorm(K, 0, 2)
  Beta[,1] = beta = betan = rep(0,P)
  Pi_beta[,1] = pi_beta = pi_betan = rep(0.1,P)
  logpi_beta = logpi_betan = log(pi_beta)
  Omega_[1:K,1:K,1] = Omega = Omegan = matrix(rgamma(K*K, 1, 1), K, K)
  omega_[1:K,1] = omega = omegan = rgamma(K, 1, 1)
  U[,,1] = u = sample(c(1:K), N*Ti, replace = TRUE, prob = rep(1/K, K))
  u = matrix(u, N, Ti)
  
  if(K==1){
    Pi_[,,1] = Pi = 1
    pi_[,1] = pi = 1
  }else{
    Pi_[1:K,1:K,1] = Pi = Omega_[1:K,1:K,1]/rowSums(Omega_[1:K,1:K,1])
    pi_[1:K,1] = pi = omega_[1:K,1]/sum(omega_[1:K,1])
  }
  
  #Hyperparameters
  mu_bt = 0
  sig_bt = 1
  
  #Hyperparameters
  mu_gm = 0
  sig_gm = 1
  
  #Spike_slab
  spike = c(1, 0.5, 0.1, 0.05, 0.01, 0.005)
  lp = length(spike)
  
  spike2 = c(1, 0.5, 0.1, 0.05, 0.01, 0.005)
  lp2 = length(spike2)
  
  sig_spike_ = rep(spike, c(rep(R/(2*(lp-1)), lp-1), R/2))
  sig_spike_2 = rep(spike2, c(rep(R/(2*(lp2-1)), lp2-1), R/2))
  
  a_pi = 1/2
  b_pi = 1/2
  
  mu_al = 0
  sig_al = 1
  hypgamma = 0.6
  vec = matrix(hypgamma, K, K)
  diag(vec) = K
  
  #Sigma RW
  sig_alp = rep(0.5, K_max)
  sig_gam = 0.5
  sig_bet = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
  sig_pi = 10
  sig_ome = 1
  sig_Ome = 0.1
  
  #Acceptance
  acc_bt = rep(0, P)
  acc_gm = 0
  acc_al = rep(0, K_max)
  acc_pi = rep(0, P)
  acc_pi_gam = 0
  acc_omega = 0
  acc_Omega = 0
  acc_com = 0
  acc_spl = 0
  acc_bir = 0
  acc_del = 0
  epsilon =  1e-10
  
  #Imputation
  Y_mcmc = array(NA, dim = c(dim(Y)[1], dim(Y)[2], R))
  
  #MCMC
  for(r in 2:R){
    
    N = as.numeric(dim(X)[1])
    Ti = as.numeric(dim(X)[3])
    for(t in 1:Ti){
      for(i in 1:N){
        if(Y_impute[i,t]==1){
          psi = alpha[u[i,t]] + beta%*%X[i,,t] 
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          Y[i,t] = rbinom(1,1,p)
        }
      }
    }
    
    Y_mcmc[,,r] = Y
    
    vec = matrix(hypgamma, K, K)
    diag(vec) = K
    
    sig_spike = sig_spike_[r]
    sig_spike2 = sig_spike_2[r]
    
    dens =  array(0, dim = c(N,Ti,K))
    for(n in 1:N){
      for(t in 1:Ti){
        for(i in 1:K){
          psi = alpha[i] + beta%*%X[n,,t] 
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          if(p == 0){
            p = p + epsilon
          }
          if(p == 1){
            p = p - epsilon
          }
          dens[n,t,i] = dbinom(Y[n,t], 1, p)
        }
      }
    }
    
    for(n in 1:N){
      for(t in Ti:1){
        if(t==Ti){
          pz <- Pi_[u[n,t-1],1:K,r-1]*dens[n,t,]
          Pz <- pz/sum(pz)
          u[n,t] = which(rmultinom(1,1,Pz)==1)
        }else if(1<t & t<Ti){
          pz <- Pi_[u[n,t-1],1:K,r-1]*dens[n,t,]*Pi_[1:K,u[n,t+1],r-1]
          Pz <- pz/sum(pz)
          u[n,t] = which(rmultinom(1,1,Pz)==1)
        }else{
          pz <- Pi_[1:K,u[n,t+1],r-1]*dens[n,t,]*pi_[1:K,r-1]
          Pz <- pz/sum(pz)
          u[n,t] = which(rmultinom(1,1,Pz)==1)
        }
      }
    }
    
    for(i in 1:K){
      alphan[i] = alpha[i] + rnorm(1, 0, sig_alp[i])
      log_ratio_alpha = log_lik(Y, u, X, alphan, beta, epsilon) + log_alpha(alphan[i], mu_al, sig_al) - log_lik(Y, u, X, alpha, beta, epsilon) - log_alpha(alpha[i], mu_al, sig_al)
      alp_alpha = min(1, exp(log_ratio_alpha))
      if(runif(1)<= alp_alpha){
        alpha[i] = alphan[i]
        acc_al[i] = acc_al[i] + 1
      }else{
        alpha[i] = alpha[i]
      }
      alphan = alpha
    }
    
    for(i in 1:P){
      betan[i] = beta[i] + rnorm(1, 0, sig_bet[i])
      log_ratio_beta = log_lik(Y, u, X, alpha, betan, epsilon) + log_beta(betan[i], mu_bt, sig_bt, sig_spike, pi_beta[i]) - log_lik(Y, u, X, alpha, beta, epsilon) - log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_beta[i])
      alp_beta = min(1, exp(log_ratio_beta))
      if(runif(1)<= alp_beta){
        beta[i] = betan[i]
        acc_bt[i] = acc_bt[i] +1
      }else{
        beta[i] = beta[i]
      }
      betan = beta
    }
    Beta[,r] = beta
    
    #Metropolis pi
    for(i in 1:P){
      if(pi_beta[i] > 1 - 1e-10){
        pi_beta[i] = 1 - 1e-10
      }else if(pi_beta[i] < 1e-10){
        pi_beta[i] = 1e-10
      }
      logpi_beta[i] = log(pi_beta[i]/(1-pi_beta[i]))
      logpi_betan[i] = logpi_beta[i] + rnorm(1, 0, sig_pi)
      pi_betan[i] = exp(logpi_betan[i])/(1+exp(logpi_betan[i])) 
      log_ratio_pi = (log_prior_pi(pi_betan[i], a_pi, b_pi) + log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_betan[i])) - (log_prior_pi(pi_beta[i], a_pi, b_pi) + log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_beta[i])) + log_jac(logpi_betan[i], logpi_beta[i])
      alp_pi = min(1,exp(log_ratio_pi))
      if(runif(1)<= alp_pi){
        pi_beta[i] = pi_betan[i]
        acc_pi[i] = acc_pi[i] +1
      }else{
        pi_beta[i] = pi_beta[i]
      }
      pi_betan = pi_beta
      logpi_betan = logpi_beta
    }
    Pi_beta[,r] = pi_beta
    
    log_omega_ = log(omega)
    log_omegan = log_omega_ + rnorm(K,0,sig_ome)
    omegan = exp(log_omegan)
    pin = omegan/sum(omegan)
    log_ratio_omega = log_latent(u, pin, Pi, K) - log_latent(u, pi, Pi, K) + log_omega(omegan, 1, 1) - log_omega(omega, 1, 1) + Jacobian_omega(omegan) - Jacobian_omega(omega)
    alp_omega = min(1, exp(log_ratio_omega))
    if(runif(1)<= alp_omega){
      omega = omegan
      pi = pin
      acc_omega = acc_omega +1
    }else{
      omega = omega
      pi = pi
    }
    
    log_Omega_ = log(Omega)
    log_Omegan = log_Omega_ + matrix(rnorm(K*K,0,sig_Ome), K, K)
    Omegan = exp(log_Omegan)
    Pin = Omegan/rowSums(Omegan)
    log_ratio_Omega = log_latent(u, pi, Pin, K) - log_latent(u, pi, Pi, K) + log_omega2(Omegan, vec) - log_omega2(Omega, vec) + Jacobian_Omega(Omegan) - Jacobian_Omega(Omega)
    alp_Omega = min(1, exp(log_ratio_Omega))
    if(runif(1)<= alp_Omega){
      Omega = Omegan
      Pi = Pin
      acc_Omega = acc_Omega +1
    }else{
      Omega = Omega
      Pi = Pi
    }
    
    #Split and Combine moves
    Omega = as.matrix(Omega)
    if(K==1){
      K = K+1
      res = 1
      pkpo = 1
      pkmo = 1/2
    }else if(K==K_max){
      K = K-1
      res = -1
      pkmo = 1
      pkpo = 1/2
    }else{
      res = sample(c(1,-1),1, 0.5)
      K = K + res
      pkpo = pkmo = 1/2
    }
    
    if(res==1){
      a = rbeta(1,10,10)
      b = rnorm(1,1/2,sqrt(0.1))
      c = runif((K-1)-1)
      d = rgamma((K-1)-1, 10, 10)
      e = rbeta(1,10,10)
      f = rgamma(1,10,10)
      g = rgamma(1,10,10)
      vec = matrix(hypgamma, K, K)
      diag(vec) = K
      j_star = sample(1:(K-1),1)
      omega_j1 = omega[j_star]*a
      omega_j2 = omega[j_star]*(1-a)
      alpha_j1 = 2*alpha[j_star]*b
      alpha_j2 = 2*alpha[j_star]*(1-b)
      if((K-1)==1){
        Omega_j1j1 = Omega[j_star, j_star]*e*f
        Omega_j1j2 = Omega[j_star, j_star]*(1-e)*f
        Omega_j2j1 = Omega[j_star, j_star]*e/f
        Omega_j2j2 = Omega[j_star, j_star]*(1-e)/f
      }else{
        Omega_hj1 = Omega[-j_star,j_star]*c
        Omega_hj2 = Omega[-j_star,j_star]*(1-c)
        Omega_j1h = Omega[j_star,-j_star]*d
        Omega_j2h = Omega[j_star,-j_star]/d
        Omega_j1j1 = Omega[j_star, j_star]*e*f
        Omega_j1j2 = Omega[j_star, j_star]*(1-e)*g
        Omega_j2j1 = Omega[j_star, j_star]*e/f
        Omega_j2j2 = Omega[j_star, j_star]*(1-e)/g
      }
      omega_split = rep(NA,K)
      omega_split[-c(j_star,j_star+1)] = omega[-j_star]
      omega_split[c(j_star,j_star+1)] = c(omega_j1,omega_j2)
      alpha_split = rep(NA,K)
      alpha_split[-c(j_star,j_star+1)] = alpha[-j_star]
      alpha_split[c(j_star,j_star+1)] = c(alpha_j1,alpha_j2)
      if((K-1)==1){
        Omega_split = matrix(c(Omega_j1j1, Omega_j2j1, Omega_j1j2, Omega_j2j2), 2, 2)
      }else{
        Omega_split = matrix(NA,K,K)
        Omega_split[-c(j_star,j_star+1),-c(j_star,j_star+1)] = Omega[-j_star,-j_star]
        Omega_split[j_star,j_star] = Omega_j1j1
        Omega_split[j_star,j_star+1] = Omega_j1j2
        Omega_split[j_star+1,j_star+1] = Omega_j2j2
        Omega_split[j_star+1,j_star] = Omega_j2j1
        Omega_split[j_star,-c(j_star,j_star+1)] = Omega_j1h 
        Omega_split[j_star+1,-c(j_star,j_star+1)] = Omega_j2h 
        Omega_split[-c(j_star,j_star+1),j_star] = Omega_hj1 
        Omega_split[-c(j_star,j_star+1),j_star+1] = Omega_hj2  
      }
      Pi_split = Omega_split/rowSums(Omega_split)
      pi_split = omega_split/sum(omega_split)
      u_split = ifelse(u <= j_star, u, u+1)
      u_split = ifelse(u_split == j_star, 0, u_split) 
      dim = sum(1*(u_split==0))
      allocation = sample(c(j_star, j_star+1), dim, replace = TRUE)
      u_split[u_split==0] = allocation
      log_P_all = log(0.5)*dim
      if(K == 2){
        Rat = log_lik(Y, u_split, X, alpha_split, beta, epsilon) - log_lik(Y, u, X, alpha, beta, epsilon) + log_latent(u_split, pi_split, Pi_split, K) - log_latent(u, pi, Pi, K-1) + log_alpha(alpha_split, mu_al, sig_al) - log_alpha(alpha, mu_al, sig_al) + log_omega(omega_split, 1, 1) - log_omega(omega, 1, 1) + log_omega2(Omega_split, vec) - log_omega2(Omega, vec[-K,-K]-diag(1,K-1))  + log(K) - sum(dgamma(c(f,g),10,10,log=TRUE)) - dnorm(b,1/2,sqrt(0.1),log=TRUE) - sum(dbeta(c(a,e),10,10,log = TRUE)) + log(pkmo) - log(pkpo) + Jacobian_split_1(omega[j_star], alpha[j_star], Omega[j_star, j_star], e, f, g) - log_P_all  
      }else{
        Rat = log_lik(Y, u_split, X, alpha_split, beta, epsilon) - log_lik(Y, u, X, alpha, beta, epsilon) + log_latent(u_split, pi_split, Pi_split, K) - log_latent(u, pi, Pi, K-1) + log_alpha(alpha_split, mu_al, sig_al) - log_alpha(alpha, mu_al, sig_al) + log_omega(omega_split, 1, 1) - log_omega(omega, 1, 1) + log_omega2(Omega_split, vec) - log_omega2(Omega, vec[-K,-K]-diag(1,K-1))  + log(K) - sum(dgamma(c(d,f,g),10,10,log=TRUE)) - dnorm(b,1/2,sqrt(0.1),log=TRUE) - sum(dbeta(c(a,e),10,10,log = TRUE)) + log(pkmo) - log(pkpo) + Jacobian_split(omega[j_star], alpha[j_star], Omega[-j_star,j_star], Omega[j_star,-j_star], Omega[j_star, j_star], d, e, f, g) - log_P_all  
      }
      alp_split = min(1,exp(Rat))
      if(runif(1)<= alp_split){
        u = u_split
        omega = omega_split
        Omega = Omega_split
        alpha = alpha_split
        Pi = Pi_split
        pi = pi_split
        K = K
        acc_spl = acc_spl +1
      }else{
        u = u
        omega = omega
        Omega = Omega
        alpha = alpha
        Pi = Pi
        pi = pi
        K = K-1
      }
    }else{
      j1 = sample(1:K+1,1)
      if(j1 == K+1){
        j1 = K
        j2 = K + 1
      }else{
        j2 = j1 +1
      }
      omega_js = omega[j1] + omega[j2]
      alpha_js = (alpha[j1]+alpha[j2])/2
      if(K == 1){
        Omega_jsjs = (Omega[j1,j1]*Omega[j2,j1])^0.5 + (Omega[j1,j2]*Omega[j2,j2])^0.5 
      }else{
        Omega_jsh = (Omega[j1,-c(j1,j2)]*Omega[j2,-c(j1,j2)])^0.5
        Omega_hjs = Omega[-c(j1,j2),j1] + Omega[-c(j1,j2),j2]
        Omega_jsjs = (Omega[j1,j1]*Omega[j2,j1])^0.5 + (Omega[j1,j2]*Omega[j2,j2])^0.5 
      }
      omega_combine = rep(NA, K)
      omega_combine[-j1] = omega[-c(j1,j2)]
      omega_combine[j1] = omega_js
      alpha_combine = rep(NA, K)
      alpha_combine[-j1] = alpha[-c(j1,j2)]
      alpha_combine[j1] = alpha_js
      if(K == 1){
        Omega_combine = Omega_jsjs
      }else{
        Omega_combine = matrix(NA,K,K)
        Omega_combine[-j1,-j1] = Omega[-c(j1,j2),-c(j1,j2)]
        Omega_combine[j1,j1] = Omega_jsjs
        Omega_combine[j1,-j1] = Omega_jsh
        Omega_combine[-j1,j1] = Omega_hjs
      }
      if(K == 1){
        Pi_combine = 1
      }else{
        Pi_combine = Omega_combine/rowSums(Omega_combine)
      }
      pi_combine = omega_combine/sum(omega_combine)
      u_combine = ifelse(u==j1 | u ==j2, j1, u)
      u_combine = ifelse(u_combine <= j1, u_combine, u_combine-1)
      
      if(K == 1){
        d = 1
      }else{
        d = sqrt((Omega[j1,-c(j1,j2)]/Omega[j2,-c(j1,j2)])^(-1))
      }
      a = omega[j1]/(omega[j1]+omega[j2])
      b = alpha[j1]/(alpha[j1]+alpha[j2])
      f = sqrt(Omega[j1,j1]/Omega[j2,j1])
      e = Omega[j1,j1]/(Omega_jsjs*f)
      g = sqrt(Omega[j1,j2]/Omega[j2,j2])
      
      dim = sum(1*(u==c(j1))) + sum(1*(u==c(j2)))
      log_P_all = log(0.5)*dim
      
      if(K == 1){
        Rat = (log_lik(Y, u, X, alpha, beta, epsilon) - log_lik(Y, u_combine, X, alpha_combine, beta, epsilon) + log_latent(u, pi, Pi, K+1) - log_latent(u_combine, pi_combine, Pi_combine, K) + log_alpha(alpha, mu_al, sig_al) - log_alpha(alpha_combine, mu_al, sig_al) + log_omega(omega, 1, 1) - log_omega(omega_combine, 1, 1) + log_omega2(Omega, vec) - log_omega2(Omega_combine, vec[-K,-K] - diag(1,K))  + log(K) - sum(dgamma(c(f,g),10,10,log=TRUE)) - dnorm(b,1/2,sqrt(0.1),log=TRUE) - sum(dbeta(c(a,e),10,10,log = TRUE)) + log(pkmo) - log(pkpo) + Jacobian_split_1(omega_js, alpha_js, Omega_jsjs, e, f, g) - log_P_all)*(-1)
      }else{
        Rat = (log_lik(Y, u, X, alpha, beta, epsilon) - log_lik(Y, u_combine, X, alpha_combine, beta, epsilon) + log_latent(u, pi, Pi, K+1) - log_latent(u_combine, pi_combine, Pi_combine, K) + log_alpha(alpha, mu_al, sig_al) - log_alpha(alpha_combine, mu_al, sig_al) + log_omega(omega, 1, 1) - log_omega(omega_combine, 1, 1) + log_omega2(Omega, vec) - log_omega2(Omega_combine, vec[-K,-K] - diag(1,K))  + log(K) - sum(dgamma(c(d,f,g),10,10,log=TRUE)) - dnorm(b,1/2,sqrt(0.1),log=TRUE) - sum(dbeta(c(a,e),10,10,log = TRUE)) + log(pkmo) - log(pkpo) + Jacobian_split(omega_js, alpha_js, Omega_hjs, Omega_jsh, Omega_jsjs, d, e, f, g) - log_P_all)*(-1)
      }
      alp_comb = min(1,exp(Rat))
      if(runif(1)<= alp_comb){
        u = u_combine
        omega = omega_combine
        Omega = Omega_combine
        alpha = alpha_combine
        Pi = Pi_combine
        pi = pi_combine
        K = K
        acc_com = acc_com +1
      }else{
        u = u
        omega = omega
        Omega = Omega
        alpha = alpha
        Pi = Pi
        pi = pi
        K = K+1
      }
    }
    
    #Birth and death moves
    if(K==1){
      K = K+1
      res = 1
      pkpo = 1
      pkmo = 1/2
    }else if(K==K_max){
      K = K-1
      res = -1
      pkmo = 1
      pkpo = 1/2
    }else{
      res = sample(c(1,-1),1, 0.5)
      K = K + res
      pkpo = pkmo = 1/2
    }
    
    if(res == 1){
      alpha_birth = c(alpha, rnorm(1, mu_al, sig_al))
      omega_birth = c(omega, rgamma(1, 1, 1))
      Omega_birth = rbind(cbind(Omega, rgamma(K-1,shape = hypgamma,1)), c(rgamma(K-1, shape = hypgamma,1), rgamma(1, shape = K, 1)))
      pi_birth = omega_birth/sum(omega_birth)
      Pi_birth = Omega_birth/rowSums(Omega_birth)
      K_0 = sum(1*table(factor(u, levels = 1:(K-1)))==0)
      Rat = log_latent(u, pi_birth, Pi_birth, K) - log_latent(u, pi, Pi, K-1) + log(K+1) - log(K_0+1)  + log(pkpo) - log(pkmo) 
      alp_bir = min(1, exp(Rat))
      if(runif(1)<= alp_bir){
        omega = omega_birth
        Omega = Omega_birth
        alpha = alpha_birth
        Pi = Pi_birth
        pi = pi_birth
        K = K
        acc_bir = acc_bir +1
      }else{
        omega = omega
        Omega = Omega
        alpha = alpha
        Pi = Pi
        pi = pi
        K = K-1
      }
      Alpha[1:K,r] = alpha
      Omega_[1:K,1:K,r] = Omega
      omega_[1:K,r] = omega
      pi_[1:K,r] = pi
      Pi_[1:K,1:K,r] = Pi
      K_state[r] = K
    }else{
      if(sum(1*table(factor(u, levels = 1:(K+1)))==0)==0){
        K = K + 1
        acc_del = acc_del
      }else{
        del = which(table(factor(u, levels = 1:(K+1)))==0)
        if(length(del) > 1){
          del = sample(c(del), 1, prob = rep(1/length(del), length(del)))
        }else{
          del = del
        }
        omega_del = omega[-del]
        Omega_del = Omega[-del,-del]
        alpha_del = alpha[-del]
        pi_del = omega_del/sum(omega_del)
        if(K == 1){
          Pi_del = 1
        }else{
          Pi_del = Omega_del/rowSums(Omega_del)
        }
        K_0 = sum(1*table(factor(u, levels = 1:(K+1)))==0)
        Rat = (log_latent(u, pi, Pi, K+1) - log_latent(u, pi_del, Pi_del, K) + log(K+1) - log(K_0+1)  + log(pkpo) - log(pkmo))*(-1) 
        alp_del = min(1, exp(Rat))
        if(runif(1)<= alp_del){
          omega = omega_del
          Omega = Omega_del
          alpha = alpha_del
          Pi = Pi_del
          pi = pi_del
          K = K
          u = u + matrix(as.vector(-1*(u>del)), N, Ti)
          acc_del = acc_del +1
        }else{
          omega = omega
          Omega = Omega
          alpha = alpha
          Pi = Pi
          pi = pi
          K = K+1
        }
      }
      Alpha[1:K,r] = alpha
      Omega_[1:K,1:K,r] = Omega
      omega_[1:K,r] = omega
      pi_[1:K,r] = pi
      Pi_[1:K,1:K,r] = Pi
      K_state[r] = K
    }
    
    #Identifiability constrint
    norm_e = Alpha[1:K,r]
    ind = order(norm_e, decreasing = TRUE)
    if(any(ind!=(1:K))){
      Alpha[1:K,r] = Alpha[ind,r] = alpha[ind]
      pi_[1:K,r] = pi_[ind,r]
      Pi_[1:K,1:K,r] = Pi_[ind,ind,r]
      omega_[1:K,r] = omega_[ind,r]
      Omega_[1:K,1:K,r] = Omega_[ind,ind,r]
      u1 = u
      old=c(1:K)
      new=ind
      u = matrix(new[match(u1,old)], N, Ti)
    }
    U[,,r] = u
    alphan = alpha
    omegan = omega
    Omegan = Omega
    Pin = Pi
    pin = pi
  }
  
  #--------------------------------
  #Bayesian HMM F-score (APP2)
  #--------------------------------
  
  prob = function(y, u, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    prob = array(NA, dim = c(N,Ti,R))
    for(r in 1:R){
      for(i in 1:N){
        for(t in 1:Ti){
          psi = alpha[u[i,t,r],r] + beta[,r]%*%x[i,,t]
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          if(p == 0){
            p = p + epsilon
          }
          if(p == 1){
            p = p - epsilon
          }
          prob[i,t,r] = p
        }
      }
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, U[,,burn_in:R], X, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = rep(NA, Tr)
  opt_thres = rep(NA, fin_it)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(Y_mcmc[,,(r+fin_it)],y_est[,,i])
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    opt_thres[r] = thres[which.max(f_score)]
  }
  
  prob_pred = function(u, x, alpha, beta, epsilon, R, N, Pi, K){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    for(r in 1:R){
      p = rep(NA, K[r])
      for(i in 1:N){
        for(v in 1:K[r]){
          psi[v] = alpha[v,r] + beta[,r]%*%x[i,] 
          p[v] = exp(psi[v])/(1+exp(psi[v]))
          if(is.nan(p[v])==TRUE){
            p[v] = 1
          }
          if(p[v] == 0){
            p[v] = p[v] + epsilon
          }
          if(p[v] == 1){
            p[v] = p[v] - epsilon
          }
        }
        P_pred[i,r] = sum(as.vector(p)*Pi[u[i,r],!is.na(Pi[u[i,r],,r]),r])  
      }
    }
    res = list(P_pred = P_pred)
  }
  
  if(ii == 29){
    y_prob_pred = prob_pred(U[,T_f,burn_in:R], X_test, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }else{
    y_prob_pred = prob_pred(U[,T_f,burn_in:R], X_test[,,1], Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }
  y_pred = matrix(NA,N,fin_it)
  for(r in 1:fin_it){ y_pred[,r] = (y_prob_pred$P_pred[,r] > opt_thres[r])*1 }
  
  Y_pred = rep(NA, N)
  
  for(i in 1:N){
    Y_pred[i] = which.max(table(factor(y_pred[i,], levels = 0:1)))[[1]]
  }
  
  Y_pred = ifelse(Y_pred==1, 0, 1)
  if(ii == 29){
    Y_test_control = Y_test
  }else{
    Y_test_control = Y_test[,1]
  }
  
  Y_pred_nona = Y_pred[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
  
  file_ = paste0("RJMCMC_spsl_and_logistic_new_imputation_50k",ii,".Rdata")
  save.image(file = file_)
  
  #--------------------------------------------
  #Classical F-score frequentist logistic model
  #--------------------------------------------
  
  Y = readRDS("Y_matrix_missing")
  X = readRDS("X_complete_new")
  
  #same scale
  X[,7,] = 100*X[,7,]
  
  #Windows
  Y_test = Y[,(T_f+ii):T_fin]
  X_test = X[,,(T_f+ii):T_fin]
  Y = Y[,ii:(T_m1+ii)]
  X = X[,,ii:(T_m1+ii)]
  
  #Classical
  mod = glm(c(Y)[which(is.na(Y)==FALSE)] ~ c(X[,1,])[which(is.na(Y)==FALSE)] + c(X[,2,])[which(is.na(Y)==FALSE)] + c(X[,3,])[which(is.na(Y)==FALSE)] + c(X[,4,])[which(is.na(Y)==FALSE)] + c(X[,5,])[which(is.na(Y)==FALSE)] + c(X[,6,])[which(is.na(Y)==FALSE)] + c(X[,7,])[which(is.na(Y)==FALSE)], family = "binomial")
  mod$coefficients[summary(mod)$coefficients[,4]>0.05] = 0
  proba = matrix(exp(mod$coefficients[1] + matrix(aperm(X, c(1, 3, 2)), ncol = P)%*%mod$coefficients[-1])/(1+exp(mod$coefficients[1] + matrix(aperm(X, c(1, 3, 2)), ncol = P)%*%mod$coefficients[-1])),n,6)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = rep(NA, Tr)
  
  for(i in 1:Tr){
    y_est[,,i] = (proba > thres[i])*1
    use = table(Y,y_est[,,i])
    precision_ = use[4]/(use[4] + use[3])
    recall_ = use[4]/(use[4] + use[2])
    f_score[i] = 2*(precision_*recall_)/(precision_ + recall_) 
  }
  opt_thres_no_bayes = thres[which.max(f_score)]
  
  if(ii == 29){
    proba_pred = matrix(exp(mod$coefficients[1] + X_test[,]%*%mod$coefficients[-1])/(1+exp(mod$coefficients[1] + X_test[,]%*%mod$coefficients[-1])),n,1)
  }else{
    proba_pred = matrix(exp(mod$coefficients[1] + X_test[,,1]%*%mod$coefficients[-1])/(1+exp(mod$coefficients[1] + X_test[,,1]%*%mod$coefficients[-1])),n,1)
  }
  y_pred_no_bayes = matrix(NA,N,fin_it)
  y_pred_no_bayes = (proba_pred > opt_thres_no_bayes)*1
  
  Y_pred_no_bayes_nona = y_pred_no_bayes[which(is.na(Y_test_control)==FALSE)]
  f_score_no_bayes = (2*table(Y_pred_no_bayes_nona,Y_test_nona)[4])/(2*table(Y_pred_no_bayes_nona,Y_test_nona)[4] + table(Y_pred_no_bayes_nona,Y_test_nona)[3] + table(Y_pred_no_bayes_nona,Y_test_nona)[2])
  
  file_ = paste0("Classical_Logistic_new_imputation_50k",ii,".Rdata")
  save.image(file = file_)
}

#---------------------------
#Bayesian HMM F-score (APP1)
#---------------------------

f_score_bayes_bart = rep(NA,29)

for(i_new in 1:29){
  load(paste0("RJMCMC_spsl_and_logistic_new_imputation_50k",i_new,".Rdata"))
  
  prob = function(y, u, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    prob = array(NA, dim = c(N,Ti,R))
    for(r in 1:R){
      for(i in 1:N){
        for(t in 1:Ti){
          psi = alpha[u[i,t,r],r] + beta[,r]%*%x[i,,t]
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          if(p == 0){
            p = p + epsilon
          }
          if(p == 1){
            p = p - epsilon
          }
          prob[i,t,r] = p
        }
      }
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, U[,,burn_in:R], X, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = matrix(NA, fin_it, Tr)
  F1b = matrix(0,fin_it,Tr)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(Y_mcmc[,,r+fin_it],y_est[,,i])
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[r,i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    F1b[r,] = f_score[r,]==max(f_score[r,],na.rm=TRUE)
  }
  
  cb1 = thres[which.max(colSums(F1b,na.rm=TRUE))]
  
  prob_pred = function(u, x, alpha, beta, epsilon, R, N, Pi, K){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    for(r in 1:R){
      p = rep(NA, K[r])
      for(i in 1:N){
        for(v in 1:K[r]){
          psi[v] = alpha[v,r] + beta[,r]%*%x[i,] 
          p[v] = exp(psi[v])/(1+exp(psi[v]))
          if(is.nan(p[v])==TRUE){
            p[v] = 1
          }
          if(p[v] == 0){
            p[v] = p[v] + epsilon
          }
          if(p[v] == 1){
            p[v] = p[v] - epsilon
          }
        }
        P_pred[i,r] = sum(as.vector(p)*Pi[u[i,r],!is.na(Pi[u[i,r],,r]),r])  
      }
    }
    res = list(P_pred = P_pred)
  }
  
  if(ii == 29){
    y_prob_pred_ = prob_pred(U[,T_f,burn_in:R], X_test, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }else{
    y_prob_pred_ = prob_pred(U[,T_f,burn_in:R], X_test[,,1], Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }
  
  y_pred_ = matrix(NA,N,fin_it)
  Y_pred_ = rep(NA, N)
  Y_pred_ = (rowMeans(y_prob_pred_$P_pred) > cb1)*1
  
  if(ii == 29){
    Y_test_control = Y_test
  }else{
    Y_test_control = Y_test[,1]
  }
  
  Y_pred_nona = Y_pred_[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes_bart[i_new] = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
  
}

save.image("Sim_RJ_Apl_new_Bart_new_imputation_50k.Rdata")

#-----------------------------------
#Bayesian HMM model averaging (APP1)
#-----------------------------------

f_score_bayes_bart_ma = rep(NA,29)

for(i_new in 1:29){
  load(paste0("RJMCMC_spsl_and_logistic_new_imputation_50k",i_new,".Rdata"))
  
  prob = function(y, u, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    prob = array(NA, dim = c(N,Ti,R))
    for(r in 1:R){
      for(i in 1:N){
        for(t in 1:Ti){
          psi = alpha[u[i,t,r],r] + beta[,r]%*%x[i,,t]
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          if(p == 0){
            p = p + epsilon
          }
          if(p == 1){
            p = p - epsilon
          }
          prob[i,t,r] = p
        }
      }
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, U[,,burn_in:R], X, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = matrix(NA, fin_it, Tr)
  F1b = matrix(0,fin_it,Tr)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(Y_mcmc[,,r+fin_it],y_est[,,i])
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[r,i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    F1b[r,] = f_score[r,]==max(f_score[r,],na.rm=TRUE)
  }
  
  cbk = rep(NA,K_max)
  F1b[is.na(F1b)]=0
  for(kk  in 1:K_max) {
    if(sum((K_state[burn_in:R]==kk)*1)==1){
      cbk[kk] = thres[which.max(F1b[K_state[burn_in:R]==kk,])]
    }else{
      cbk[kk] = thres[which.max(colSums(F1b[K_state[burn_in:R]==kk,],na.rm=TRUE))]
    }
  }
  cb1 = sum(cbk*table(factor(K_state[burn_in:R], levels = 1:K_max))/fin_it)
  
  
  prob_pred = function(u, x, alpha, beta, epsilon, R, N, Pi, K){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    for(r in 1:R){
      p = rep(NA, K[r])
      for(i in 1:N){
        for(v in 1:K[r]){
          psi[v] = alpha[v,r] + beta[,r]%*%x[i,] 
          p[v] = exp(psi[v])/(1+exp(psi[v]))
          if(is.nan(p[v])==TRUE){
            p[v] = 1
          }
          if(p[v] == 0){
            p[v] = p[v] + epsilon
          }
          if(p[v] == 1){
            p[v] = p[v] - epsilon
          }
        }
        P_pred[i,r] = sum(as.vector(p)*Pi[u[i,r],!is.na(Pi[u[i,r],,r]),r])  
      }
    }
    res = list(P_pred = P_pred)
  }
  
  if(ii == 29){
    y_prob_pred_ = prob_pred(U[,T_f,burn_in:R], X_test, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }else{
    y_prob_pred_ = prob_pred(U[,T_f,burn_in:R], X_test[,,1], Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }
  
  
  y_pred_ = matrix(NA,N,fin_it)
  Y_pred_ = rep(NA, N)
  y_prob_nma = matrix(NA,K_max,N)
  for(kk  in 1:K_max){
    if(sum((K_state[burn_in:R]==kk)*1)==1){
      y_prob_nma[kk,] = y_prob_pred_$P_pred[,K_state[burn_in:R]==kk]
    }else{
      y_prob_nma[kk,] = rowMeans(y_prob_pred_$P_pred[,K_state[burn_in:R]==kk])
    }
  }
  y_prob_nma[is.na(y_prob_nma)] = 0
  y_prob_ma = t(y_prob_nma)%*%table(factor(K_state[burn_in:R], levels = 1:K_max))/fin_it
  Y_pred_ = (y_prob_ma > cb1)*1
  
  
  if(ii == 29){
    Y_test_control = Y_test
  }else{
    Y_test_control = Y_test[,1]
  }
  
  Y_pred_nona = Y_pred_[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes_bart_ma[i_new] = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
  
}

save.image("Sim_RJ_Apl_new_Bart_new_imputation_50k_model_average.Rdata")

#--------------------------
#Bayesian HMM model averaging (APP2)
#--------------------------

f_score_bayes_tanc_ma = rep(NA,29)

for(i_new in 1:29){
  load(paste0("RJMCMC_spsl_and_logistic_new_imputation_50k",i_new,".Rdata"))
  
  prob = function(y, u, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    prob = array(NA, dim = c(N,Ti,R))
    for(r in 1:R){
      for(i in 1:N){
        for(t in 1:Ti){
          psi = alpha[u[i,t,r],r] + beta[,r]%*%x[i,,t]
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          if(p == 0){
            p = p + epsilon
          }
          if(p == 1){
            p = p - epsilon
          }
          prob[i,t,r] = p
        }
      }
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, U[,,burn_in:R], X, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = rep(NA, Tr)
  opt_thres = rep(NA, fin_it)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(Y_mcmc[,,(r+fin_it)],y_est[,,i])
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    opt_thres[r] = thres[which.max(f_score)]
  }
  
  prob_pred = function(u, x, alpha, beta, epsilon, R, N, Pi, K){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    for(r in 1:R){
      p = rep(NA, K[r])
      for(i in 1:N){
        for(v in 1:K[r]){
          psi[v] = alpha[v,r] + beta[,r]%*%x[i,] 
          p[v] = exp(psi[v])/(1+exp(psi[v]))
          if(is.nan(p[v])==TRUE){
            p[v] = 1
          }
          if(p[v] == 0){
            p[v] = p[v] + epsilon
          }
          if(p[v] == 1){
            p[v] = p[v] - epsilon
          }
        }
        P_pred[i,r] = sum(as.vector(p)*Pi[u[i,r],!is.na(Pi[u[i,r],,r]),r])  
      }
    }
    res = list(P_pred = P_pred)
  }
  
  if(ii == 29){
    y_prob_pred = prob_pred(U[,T_f,burn_in:R], X_test, Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }else{
    y_prob_pred = prob_pred(U[,T_f,burn_in:R], X_test[,,1], Alpha[,burn_in:R], Beta[,burn_in:R], epsilon, fin_it, N, Pi_[,,burn_in:R], K_state[burn_in:R])
  }
  y_pred = matrix(NA,N,fin_it)
  for(r in 1:fin_it){ y_pred[,r] = (y_prob_pred$P_pred[,r] > opt_thres[r])*1 }
  
  Y_pred = rep(NA, N)
  y_pred_k = matrix(NA,N,K_max)
  for(kk in 1:K_max) {
    if(sum((K_state[burn_in:R]!=kk)*1)==fin_it){
      y_pred_k[,kk] = 0
    }else if(sum((K_state[burn_in:R]==kk)*1)==1){
      y_pred_k[,kk] = y_pred[,K_state[burn_in:R]==kk]/sum((K_state[burn_in:R]==kk)*1)
    }else{
      y_pred_k[,kk] = rowSums(y_pred[,K_state[burn_in:R]==kk])/sum((K_state[burn_in:R]==kk)*1)}
  }
  y_pred_ma = y_pred_k%*%(table(factor(K_state[burn_in:R], levels = 1:K_max))/fin_it)
  
  Y_pred = ifelse(y_pred_ma>0.5, 1, 0)
  
  if(ii == 29){
    Y_test_control = Y_test
  }else{
    Y_test_control = Y_test[,1]
  }
  
  Y_pred_nona = Y_pred[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes_tanc_ma[i_new] = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
}

save.image("Sim_RJ_Apl_new_Tanc_new_imputation_50k_model_average.Rdata")

#--------------------
#Bayesian logistic
#--------------------

#Set WD
#setwd("~/BayesianF_score")

Y = readRDS("Y_matrix_missing")
X = readRDS("X_complete_new")

#Library
library("LaplacesDemon")
library("abind")

#Functions
log_prior_pi = function(x, a, b){
  dbeta(x, a, b, log = TRUE)
}

log_jac = function(phin, phi){
  log(exp(phin)/(1+exp(phin))^2) - log(exp(phi)/(1+exp(phi))^2)
}

log_beta = function(beta, mu, sigma, sigma2, pi_beta_){
  log(pi_beta_*dlaplace(beta, mu, sigma) + (1-pi_beta_)*dlaplace(beta, mu, sigma2))
}

log_alpha = function(alpha, mu, sigma){
  sum(dnorm(alpha, mu, sigma, log = TRUE))
}

log_lik = function(y, x, alpha, beta, epsilon){
  N = as.numeric(dim(x)[1])
  Ti = as.numeric(dim(x)[3])
  P = as.numeric(dim(X)[2])
  x = matrix(aperm(x, c(1, 3, 2)), ncol = P)
  psi = alpha + x%*%beta
  p = exp(psi)/(1+exp(psi))
  if(sum(is.nan(p))>0){
    p[is.nan(p)] = 1
  }
  if(sum(p == 0)>0){
    p[p == 0] = epsilon
  }
  if(sum(p == 1)>0){
    p[p == 1] = 1 - epsilon
  }
  sum(dbinom(c(y), 1, p, log = TRUE))
}

T_fin = 35
T_m1 = 5
T_f = T_m1+1

for(ii in 22:29){
  Y = readRDS("Y_matrix_missing")
  X = readRDS("X_complete_new")
  
  #Same scale
  X[,7,] = 100*X[,7,]
  
  Y_test = Y[,(T_f+ii):T_fin]
  X_test = X[,,(T_f+ii):T_fin]
  Y = Y[,ii:(T_m1+ii)]
  X = X[,,ii:(T_m1+ii)]
  n = N = as.numeric(nrow(Y))
  Time = Ti = as.numeric(ncol(Y))
  
  Y_impute = matrix(0,N,Ti)
  for(t in 1:(T_m1+1)){
    Y_impute[which(is.na(Y[,t])),t] = 1  
  }
  
  #Number of predictors
  P = as.numeric(dim(X)[2])
  
  #Iterations
  R = 50000 
  
  #Storage
  Alpha = rep(NA, R)
  Beta = matrix(NA, P, R)
  Pi_beta = matrix(NA, P, R) 
  
  #Initial values
  Alpha[1] = alpha = alphan = rnorm(1, 0, 1)
  Beta[,1] = beta = betan = rep(0,P)
  Pi_beta[,1] = pi_beta = pi_betan = rep(0.1,P)
  logpi_beta = logpi_betan = log(pi_beta)
  
  #Hyperparameters
  mu_bt = 0
  sig_bt = 1
  
  #Spike_slab
  spike = c(1, 0.5, 0.1, 0.05, 0.01, 0.005)
  lp = length(spike)
  
  spike2 = c(1, 0.5, 0.1, 0.05, 0.01, 0.005)
  lp2 = length(spike2)
  
  sig_spike_ = rep(spike, c(rep(R/(2*(lp-1)), lp-1), R/2))
  sig_spike_2 = rep(spike2, c(rep(R/(2*(lp2-1)), lp2-1), R/2))
  
  a_pi = 1/2
  b_pi = 1/2
  
  mu_al = 0
  sig_al = 1
  
  #Sigma RW
  sig_alp = c(0.2)
  sig_bet = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2)
  sig_pi = 10
  
  #Acceptance
  acc_bt = rep(0, P)
  acc_al = c(0)
  acc_pi = rep(0, P)
  epsilon =  1e-10
  
  Y_mcmc = array(NA, dim = c(dim(Y)[1], dim(Y)[2], R))
  
  #MCMC algorithm
  for(r in 2:R){
    
    N = as.numeric(dim(X)[1])
    Ti = as.numeric(dim(X)[3])
    sig_spike = sig_spike_[r]
    sig_spike2 = sig_spike_2[r]
    
    for(t in 1:Ti){
      for(i in 1:N){
        if(Y_impute[i,t]==1){
          psi = alpha + beta%*%X[i,,t] 
          p = exp(psi)/(1+exp(psi))
          if(is.nan(p)==TRUE){
            p = 1
          }
          Y[i,t] = rbinom(1,1,p)
        }
      }
    }
    
    Y_mcmc[,,r] = Y
    
    #Metropolis alpha
    alphan = alpha + rnorm(1, 0, sig_alp)
    log_ratio_alpha = log_lik(Y, X, alphan, beta, epsilon) + log_alpha(alphan, mu_al, sig_al) - log_lik(Y, X, alpha, beta, epsilon) - log_alpha(alpha, mu_al, sig_al)
    alp_alpha = min(1, exp(log_ratio_alpha))
    if(runif(1)<= alp_alpha){
      alpha = alphan
      acc_al = acc_al + 1
    }else{
      alpha = alpha
    }
    Alpha[r] = alpha
    
    #Metropolis beta
    for(i in 1:P){
      betan[i] = beta[i] + rnorm(1, 0, sig_bet[i])
      log_ratio_beta = log_lik(Y, X, alpha, betan, epsilon) + log_beta(betan[i], mu_bt, sig_bt, sig_spike, pi_beta[i]) - log_lik(Y, X, alpha, beta, epsilon) - log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_beta[i])
      alp_beta = min(1, exp(log_ratio_beta))
      if(runif(1)<= alp_beta){
        beta[i] = betan[i]
        acc_bt[i] = acc_bt[i] +1
      }else{
        beta[i] = beta[i]
      }
      betan = beta
    }
    Beta[,r] = beta
    
    #Metropolis pi
    for(i in 1:P){
      if(pi_beta[i] > 1 - 1e-10){
        pi_beta[i] = 1 - 1e-10
      }else if(pi_beta[i] < 1e-10){
        pi_beta[i] = 1e-10
      }
      logpi_beta[i] = log(pi_beta[i]/(1-pi_beta[i]))
      logpi_betan[i] = logpi_beta[i] + rnorm(1, 0, sig_pi)
      pi_betan[i] = exp(logpi_betan[i])/(1+exp(logpi_betan[i])) 
      log_ratio_pi = (log_prior_pi(pi_betan[i], a_pi, b_pi) + log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_betan[i])) - (log_prior_pi(pi_beta[i], a_pi, b_pi) + log_beta(beta[i], mu_bt, sig_bt, sig_spike, pi_beta[i])) + log_jac(logpi_betan[i], logpi_beta[i])
      alp_pi = min(1,exp(log_ratio_pi))
      if(runif(1)<= alp_pi){
        pi_beta[i] = pi_betan[i]
        acc_pi[i] = acc_pi[i] +1
      }else{
        pi_beta[i] = pi_beta[i]
      }
      pi_betan = pi_beta
      logpi_betan = logpi_beta
    }
    Pi_beta[,r] = pi_beta
    
  }
  
  #--------------------------------
  #Bayesian logistic F-score (APP2)
  #--------------------------------
  
  prob = function(y, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    P = as.numeric(dim(X)[2])
    prob = array(NA, dim = c(N,Ti,R))
    x = matrix(aperm(x, c(1, 3, 2)), ncol = P)
    for(r in 1:R){
      psi = alpha[r] + x%*%beta[,r]
      p = exp(psi)/(1+exp(psi))
      if(sum(is.nan(p))>0){
        p[is.nan(p)] = 1
      }
      if(sum(p == 0)>0){
        p[p == 0] = epsilon
      }
      if(sum(p == 1)>0){
        p[p == 1] = 1 - epsilon
      }
      prob[,,r] = matrix(p,N,Ti)
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, X, Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = rep(NA, Tr)
  opt_thres = rep(NA, fin_it)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(factor(Y_mcmc[,,(r+fin_it)], levels= 0:1), factor(y_est[,,i], levels= 0:1))
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    opt_thres[r] = thres[which.max(f_score)]
  }
  
  prob_pred = function(x, alpha, beta, epsilon, R, N){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    x = matrix(aperm(x, c(1, 2)), ncol = P)
    for(r in 1:R){
      psi = alpha[r] + x%*%beta[,r]
      p = exp(psi)/(1+exp(psi))
      if(sum(is.nan(p))>0){
        p[is.nan(p)] = 1
      }
      if(sum(p == 0)>0){
        p[p == 0] = epsilon
      }
      if(sum(p == 1)>0){
        p[p == 1] = 1 - epsilon
      }
      prob[,r] = matrix(p,N,1)
    }
    P_pred = prob
    res = list(P_pred = P_pred)
  }
  
  if(ii==29){
    y_prob_pred = prob_pred(X_test, Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it, n)
  }else{
    y_prob_pred = prob_pred(X_test[,,1], Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it, n)
  }
  y_pred = matrix(NA,N,fin_it)
  for(r in 1:fin_it){ y_pred[,r] = (y_prob_pred$P_pred[,r] > opt_thres[r])*1 }
  Y_pred = rep(NA, N)
  for(i in 1:N){
    Y_pred[i] = which.max(table(factor(y_pred[i,], levels = 0:1)))[[1]]
  }
  Y_pred = ifelse(Y_pred==1, 0, 1)
  Y_test_control = Y_test
  
  Y_pred_nona = Y_pred[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
  
  file_ = paste0("Bayesian_logistic_pig_new_imputation_50k",ii,".Rdata")
  save.image(file = file_)
}

#--------------------------------
#Bayesian logistic F-score (APP1)
#--------------------------------

f_score_bayes_bart = rep(NA,29)

for(i_new in 1:29){
  load(paste0("Bayesian_logistic_pig_new_imputation_50k",i_new,".Rdata"))
  
  prob = function(y, x, alpha, beta, epsilon, R){
    N = as.numeric(nrow(y))
    Ti = as.numeric(ncol(y))
    P = as.numeric(dim(X)[2])
    prob = array(NA, dim = c(N,Ti,R))
    x = matrix(aperm(x, c(1, 3, 2)), ncol = P)
    for(r in 1:R){
      psi = alpha[r] + x%*%beta[,r]
      p = exp(psi)/(1+exp(psi))
      if(sum(is.nan(p))>0){
        p[is.nan(p)] = 1
      }
      if(sum(p == 0)>0){
        p[p == 0] = epsilon
      }
      if(sum(p == 1)>0){
        p[p == 1] = 1 - epsilon
      }
      prob[,,r] = matrix(p,N,Ti)
    }
    res = list(prob = prob)
  }
  
  fin_it = 25000
  burn_in = 25001
  
  probabilit = prob(Y, X, Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it)
  
  thres = seq(0.01,0.99,0.01)
  Tr = length(thres)
  
  y_est = array(NA, dim = c(N, Ti, Tr))
  f_score = matrix(NA, fin_it, Tr)
  F1b = matrix(0,fin_it,Tr)
  
  for(r in 1:fin_it){
    for(i in 1:Tr){
      y_est[,,i] = (probabilit$prob[,,r] > thres[i])*1
      use = table(Y_mcmc[,,(r+fin_it)],y_est[,,i])
      precision_ = use[4]/(use[4] + use[3])
      recall_ = use[4]/(use[4] + use[2])
      f_score[r,i] = 2*(precision_*recall_)/(precision_ + recall_) 
    }
    F1b[r,] = f_score[r,]==max(f_score[r,],na.rm=TRUE)
  }
  
  cb1 = thres[which.max(colSums(F1b,na.rm=TRUE))]
  
  prob_pred = function(x, alpha, beta, epsilon, R, N){
    prob = matrix(NA, N, R)
    P_pred = matrix(NA, N, R)
    x = matrix(aperm(x, c(1, 2)), ncol = P)
    for(r in 1:R){
      psi = alpha[r] + x%*%beta[,r]
      p = exp(psi)/(1+exp(psi))
      if(sum(is.nan(p))>0){
        p[is.nan(p)] = 1
      }
      if(sum(p == 0)>0){
        p[p == 0] = epsilon
      }
      if(sum(p == 1)>0){
        p[p == 1] = 1 - epsilon
      }
      prob[,r] = matrix(p,N,1)
    }
    P_pred = prob
    res = list(P_pred = P_pred)
  }
  
  if(ii==29){
    y_prob_pred_ = prob_pred(X_test, Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it, n)
  }else{
    y_prob_pred_ = prob_pred(X_test[,,1], Alpha[burn_in:R], Beta[,burn_in:R], epsilon, fin_it, n)
  }
  
  y_pred_ = matrix(NA,N,fin_it)
  Y_pred_ = rep(NA, N)
  Y_pred_ = (rowMeans(y_prob_pred_$P_pred) > cb1)*1
  
  if(ii == 29){
    Y_test_control = Y_test
  }else{
    Y_test_control = Y_test[,1]
  }
  
  Y_pred_nona = Y_pred_[which(is.na(Y_test_control)==FALSE)]
  Y_test_nona = Y_test_control[which(is.na(Y_test_control)==FALSE)]
  
  f_score_bayes_bart[i_new] = (2*table(Y_pred_nona,Y_test_nona)[4])/(2*table(Y_pred_nona,Y_test_nona)[4] + table(Y_pred_nona,Y_test_nona)[3] + table(Y_pred_nona,Y_test_nona)[2])
  
}

save.image("Sim_LogisticBayes_Apl_new_Bart_new_imputation_50k.Rdata")

#--------------------
#Results
#--------------------

#Set WD with results
#setwd("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score/Res_sim/Pigini_Data_new_imputation_including_Yr_50k")

#F-score values
F_score_bayes_2 = rep(NA, 29)
F_score_no_bayes = rep(NA, 29)
F_score_bayes_log2 = rep(NA, 29)
F_score_bayes = rep(NA, 29)
F_score_bayes_log = rep(NA, 29)

for(ii in 1 :29){
  load(paste0("Classical_Logistic_new_imputation_50k",ii,".Rdata"))
  F_score_no_bayes[ii] = f_score_no_bayes
}

for(ii in 1 :29){
  load(paste0("RJMCMC_spsl_and_logistic_new_imputation_50k",ii,".Rdata"))
  F_score_bayes_2[ii] = f_score_bayes
}

#Set NA as 0
F_score_bayes_2[is.na(F_score_bayes_2)] = 0
F_score_no_bayes[is.na(F_score_no_bayes)] = 0

for(ii in 1 :29){
  load(paste0("Bayesian_logistic_pig_new_imputation_50k",ii,".Rdata"))
  F_score_bayes_log2[ii] = f_score_bayes
}

F_score_bayes_log2[is.na(F_score_bayes_log2)] = 0

#Load the single Rdata
load("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score/Res_sim/Pigini_Data_new_imputation_including_Yr_50k/Sim_RJ_Apl_new_Bart_new_imputation_50k.Rdata")

F_score_bayes = f_score_bayes_bart
F_score_bayes[is.na(F_score_bayes)] = 0

load("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score/Res_sim/Pigini_Data_new_imputation_including_Yr_50k/Sim_LogisticBayes_Apl_new_Bart_new_imputation_50k.Rdata")

F_score_bayes_log = f_score_bayes_bart
F_score_bayes_log[is.na(F_score_bayes_log)] = 0

#Final values multimodel
F_score_bayes #APP1 Bayesian HMM multimodel
F_score_bayes_2 #APP2 Bayesian HMM multimodel
F_score_bayes_log #APP1 Bayesian logistic
F_score_bayes_log2 #APP2 Bayesian logistic
F_score_no_bayes #Classical logistic

#F-score values
F_score_bayes_2 = rep(NA, 29)
F_score_no_bayes = rep(NA, 29)
F_score_bayes_log2 = rep(NA, 29)
F_score_bayes = rep(NA, 29)
F_score_bayes_log = rep(NA, 29)

for(ii in 1 :29){
  load(paste0("Classical_Logistic_new_imputation_50k",ii,".Rdata"))
  F_score_no_bayes[ii] = f_score_no_bayes
}


load(paste0("Sim_RJ_Apl_new_Tanc_new_imputation_50k_model_average.Rdata"))
F_score_bayes_2 = f_score_bayes_tanc_ma

#Set NA as 0
F_score_bayes_2[is.na(F_score_bayes_2)] = 0
F_score_no_bayes[is.na(F_score_no_bayes)] = 0

for(ii in 1 :29){
  load(paste0("Bayesian_logistic_pig_new_imputation_50k",ii,".Rdata"))
  F_score_bayes_log2[ii] = f_score_bayes
}

F_score_bayes_log2[is.na(F_score_bayes_log2)] = 0

#Load the single Rdata
load("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score/Res_sim/Pigini_Data_new_imputation_including_Yr_50k/Sim_RJ_Apl_new_Bart_new_imputation_50k_model_average.Rdata")

F_score_bayes = f_score_bayes_bart_ma
F_score_bayes[is.na(F_score_bayes)] = 0

load("~/Desktop/Phd_proj/Repository lavori/Bayesian_F_score/Res_sim/Pigini_Data_new_imputation_including_Yr_50k/Sim_LogisticBayes_Apl_new_Bart_new_imputation_50k.Rdata")

F_score_bayes_log = f_score_bayes_bart
F_score_bayes_log[is.na(F_score_bayes_log)] = 0

#Final values 
F_score_bayes #APP1 Bayesian HMM multimodel
F_score_bayes_2 #APP2 Bayesian HMM multimodel
F_score_bayes_log #APP1 Bayesian logistic
F_score_bayes_log2 #APP2 Bayesian logistic
F_score_no_bayes #Classical logistic






