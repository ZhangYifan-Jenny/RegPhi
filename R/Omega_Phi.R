################## update Omega ##############################################################
Gibbs_Omega <- function(y, mu = 0, la, phi, inv.psx, N, K) {
    # y=Y;ly=LD;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + inv.phi
    SIG <- chol2inv(chol(ISG))
    Ycen <- y - mu
    mean <- SIG %*% t(la) %*% inv.psx %*% Ycen

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

######## update PHI ########################################################
MH_PHI <- function(phi, ome, N, K, s0,std) {
    # ome=Omega;ph0=PHI
    # v0 <- prior$v_PHI
    # s0 <- prior$s_PHI
    inv.cov <- rwish1(K + 2 + N, solve(tcrossprod(ome) + s0))
    cov <- chol2inv(chol(inv.cov))
    if(std){
      tmp <- chol2inv(chol(sqrt(diag(diag(cov)))))
      phi1 <- tmp %*% cov %*% tmp
      acc <- exp((K + 1)/2 * (log(det(phi1)) - log(det(phi))))
      acid <- acc > runif(1)
      phi <- phi1 * acid + phi * (1 - acid)
    }else{phi<-cov}
    
    # inv.CPH<-chol2inv(chol(CPH))
    return(phi)
}
######## end of update PHI #################################################

MH_PHI_LASSO <- function(phi, ome,inv.phi,const, prior) {
  # ome=Omega;ph0=PHI
  # v0 <- prior$v_PHI
  # s0 <- prior$s_PHI

  J <- const$J
  K <- const$K
  N <- const$N

  upperind_phi <- const$upperind_phi
  lowerind_phi <- const$lowerind_phi
  ind_nod_phi <- const$ind_nod_phi

  a_gams <- prior$a_gams
  b_gams <- prior$b_gams
  S <- ome %*% t(ome)  #k*k
  #temp <- y - mu - la %*% ome  # J*N
  #S <- temp %*% t(temp)  # J*J
  
  # sample gammas
  a_gams <- 1
  b_gams <- 0.1
  apost <- a_gams + K * (K + 1)/2
  bpost <- b_gams + sum(abs(inv.phi))/2  # C is the presicion matrix
  gammas <- rgamma(1, shape = apost, rate = bpost)
  
  # sample tau off-diagonal
  Cadj <- pmax(abs(inv.phi[upperind_phi]), 10^(-6))
  mu_p <- pmin(gammas/Cadj, 10^12)
  gammas_p <- gammas^2
  len <- length(mu_p)
  taus_tmp <- 1/rinvgauss1(len, mean = mu_p, dispersion = 1/gammas_p)
  
  taus <- matrix(0, K, K)
  taus[upperind_phi] <- taus[lowerind_phi] <- taus_tmp
  
  # sample PSX and inv(PSX)
  for (i in 1:K) {
    ind_noi <- ind_nod_phi[, i]
    Sig11 <- phi[ind_noi, ind_noi]
    Sig12 <- phi[ind_noi, i]
    invC11 <- Sig11 - Sig12 %*% t(Sig12)/phi[i, i]
    Ci <- (S[i, i] + gammas) * invC11 + diag(1/taus[ind_noi, i])
    Sigma <- chol2inv(chol(Ci))
    mu_i <- -Sigma %*% S[ind_noi, i]
    beta <- mvrnorm(1, mu_i, Sigma)
    inv.phi[ind_noi, i] <- beta
    inv.phi[i, ind_noi] <- beta
    gam <- rgamma(1, shape = N/2 + 1, rate = (S[i, i] + gammas)/2)
    inv.phi[i, i] <- gam + t(beta) %*% invC11 %*% beta
    
    # below updating covariance matrix according to one-column change of precision matrix
    invC11beta <- invC11 %*% beta
    phi[ind_noi, ind_noi] <- invC11 + invC11beta %*% t(invC11beta)/gam
    Sig12 <- -invC11beta/gam
    phi[ind_noi, i] <- Sig12
    phi[i, ind_noi] <- t(Sig12)
    phi[i, i] <- 1/gam
    
  }  # end of i

  return(list(phi = phi, inv.phi = inv.phi, gammas = gammas))
}

MH_PHI_horse <- function(phi, ome,inv.phi,const,Lambda_sq ,Nu,xi,regdiag=T) {

  J <- const$J
  K <- const$K
  N <- const$N

  ind_nod_phi <- const$ind_nod_phi

  S <- ome %*% t(ome)  #k*k
  #temp <- y - mu - la %*% ome  # J*N
  #S <- temp %*% t(temp)  # J*J
  
  ### sample tau_sq and xi
  omega_vector = inv.phi[lower.tri(inv.phi)]
  lambda_sq_vector = Lambda_sq[lower.tri(Lambda_sq)]
  rate = 1/xi + sum(omega_vector^2/(2*lambda_sq_vector));
  tau_sq = rgamma(1, shape =  K * (K + 1)/2, rate =  rate)  ## inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
  xi = rgamma(1, shape = 1, rate =  1/(1+1/tau_sq))    ## inv gamma w/ shape=1, rate=1+1/tau_sq
  if(regdiag){lamdiag<-tau_sq}else{lamdiag<-0}
  # sample PHI and inv(PHI)
  for (i in 1:K) {
    ind_noi <- ind_nod_phi[, i]
    # tau_temp1<-tau[ind_noi,i]
    Sig11 <- phi[ind_noi, ind_noi]
    Sig12 <- phi[ind_noi, i]
    
    #### sample beta and gamma
    invC11 <- Sig11 - Sig12 %*% t(Sig12)/phi[i, i]
    Ci <- (S[i, i]+lamdiag) * invC11 + diag(1/(Lambda_sq[ind_noi, i]*tau_sq))   ##lasso is diag(tau)
    Sigma <- chol2inv(chol(Ci))
    mu_i <- -Sigma %*% S[ind_noi, i]
    beta <- mvrnorm(1, mu_i, Sigma)
    
    ### upadate inverse covariance matrix
    inv.phi[ind_noi, i] <- beta
    inv.phi[i, ind_noi] <- beta
    gam <- rgamma(1, shape = N/2 + 1, rate = (S[i, i]+lamdiag)/2)
    inv.phi[i, i] <- gam + t(beta) %*% invC11 %*% beta
    
    # below updating covariance matrix according to one-column change of precision matrix
    invC11beta <- invC11 %*% beta
    phi[ind_noi, ind_noi] <- invC11 + invC11beta %*% t(invC11beta)/gam
    Sig12 <- -invC11beta/gam
    phi[ind_noi, i] <- Sig12
    phi[i, ind_noi] <- t(Sig12)
    phi[i, i] <- 1/gam
    
    ###sample lambda and nu
    rate = t(beta) %*% beta /(2*tau_sq)+1./Nu[ind_noi, i];
    Lambda_sq[i,ind_noi]=Lambda_sq[ind_noi, i] = rgamma(1, shape = 1, rate = 1/rate);    # random inv gamma with shape=1, rate=rate
    Nu[i,ind_noi]=Nu[ind_noi, i] = rgamma( 1,shape = 1, rate = 1./(1+1/Lambda_sq[ind_noi, i]));    # random inv gamma with shape=1, rate=1+1/lambda_sq_12
  }  
  return(list(phi = phi, inv.phi = inv.phi, Lambda_sq = Lambda_sq,Nu=Nu,xi=xi))
}

MH_PHI_ssp <- function(phi, ome,inv.phi,const,pii,taus,lambda=1,V0,V1) {
  
  J <- const$J
  K <- const$K
  N <- const$N
  
  ind_nod_phi <- const$ind_nod_phi
  
  S <- ome %*% t(ome)  #k*k
  #temp <- y - mu - la %*% ome  # J*N
  #S <- temp %*% t(temp)  # J*J
  # sample PSX and inv(PSX)
  for (i in 1:K) {
    
    ind_noi <- ind_nod_phi[, i]
    # tau_temp1<-tau[ind_noi,i]
    Sig11 <- phi[ind_noi, ind_noi]
    Sig12 <- phi[ind_noi, i]
    invC11 <- Sig11 - Sig12 %*% t(Sig12)/phi[i, i]
    # Ci<-(S[i,i]+gammas)*invC11+diag(1/tau_temp1)
    Ci <- (S[i, i] + lambda) * invC11 + diag(1/taus[ind_noi, i])
    Sigma <- chol2inv(chol(Ci))
    mu_i <- -Sigma %*% S[ind_noi, i]
    beta <- mvrnorm(1, mu_i, Sigma)
    inv.phi[ind_noi, i] <- beta
    inv.phi[i, ind_noi] <- beta
    gam <- rgamma(1, shape = N/2 + 1, rate = (S[i, i] + lambda)/2)
    inv.phi[i, i] <- gam + t(beta) %*% invC11 %*% beta
    
    # below updating covariance matrix according to one-column change of precision matrix
    invC11beta <- invC11 %*% beta
    phi[ind_noi, ind_noi] <- invC11 + invC11beta %*% t(invC11beta)/gam
    Sig12 <- -invC11beta/gam
    phi[ind_noi, i] <- Sig12
    phi[i, ind_noi] <- t(Sig12)
    phi[i, i] <- 1/gam
    
    # Extracting v0 and v1
    v0 <- V0[ind_noi, i]
    v1 <- V1[ind_noi, i]
    w1 <- -0.5 * log(v0) - 0.5 * beta^2 / v0 + log(1 - pii)
    w2 <- -0.5 * log(v1) - 0.5 * beta^2 / v1 + log(pii)
    w_max <- apply(cbind(w1, w2), 1, max)
    w <- exp(w2 - w_max) / rowSums(exp(cbind(w1, w2) - matrix(rep(w_max, 2), nrow = length(w_max), byrow = TRUE)))
    z <- as.numeric(runif(K - 1) < w)
    v <- v0
    v[z == 1] <- v1[z == 1]
    taus[ind_noi, i] <- v
    taus[i, ind_noi] <- v

  }  # end of i, sample Sig and C=inv(Sig)

  return(list(phi = phi, inv.phi = inv.phi, taus = taus))
}

######## update PHI ########################################################
# old version; all prior needed
MH_PHI1 <- function(phi, ome, N, K, prior) {
    # ome=Omega;ph0=PHI
    v0 <- prior$v_PHI
    s0 <- prior$s_PHI
    S<-solve(tcrossprod(ome) + s0)
    inv.cov <- rwish1(v0 + N, S)
    cov <- chol2inv(chol(inv.cov))

    tmp <- chol2inv(chol(sqrt(diag(diag(cov)))))
    phi1 <- tmp %*% cov %*% tmp
    acc <- exp((K + 1)/2 * (log(det(phi1)) - log(det(phi))))
    acid <- acc > runif(1)
    phi <- phi1 * acid + phi * (1 - acid)
    # inv.CPH<-chol2inv(chol(CPH))
    return(phi)
}
######## end of update PHI #################################################

################## update Omega ##############################################################
# allow factor mean to be specified (i.e., no zero)
Gibbs_Omega1 <- function(y, mu = 0, la, phi, inv.psx, N, K,int) {
    # y=Y;la=LA;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    # inv.phi <- solve(phi)
    inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + inv.phi
    SIG <- chol2inv(chol(ISG))
    # SIG <- solve(ISG)
    Ycen <- y - mu
    mean <- SIG %*% (t(la) %*% inv.psx %*% Ycen + inv.phi %*% int)

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

################## update Omega ##############################################################
# allow factor mean and variance to be specified
Gibbs_Omega2 <- function(y, mu = 0, la, inv.psx, N, K,m.add,s.add) {
    # y=Y;la=LA;psx=PSX;mu=0;inv.psx=inv.PSX
    inv.sqrt.psx <- chol(inv.psx)
    # inv.phi <- solve(phi)
    # inv.phi <- chol2inv(chol(phi))
    ISG <- crossprod(inv.sqrt.psx %*% la) + s.add
    SIG <- chol2inv(chol(ISG))
    # SIG <- solve(ISG)
    Ycen <- y - mu
    mean <- SIG %*% (t(la) %*% inv.psx %*% Ycen + m.add)

    CSIG <- chol(SIG)
    ome <- matrix(rnorm(N * K), K, N)
    ome <- t(t(ome) %*% CSIG) + mean

    return(ome)
}
################## end of update Omega #######################################################

