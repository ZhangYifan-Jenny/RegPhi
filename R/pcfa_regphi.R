pcfa_regphi <- function(dat, Q, LD = TRUE,cati = NULL,cand_thd = 0.2, PPMC = FALSE, burn = 5000, iter = 5000,
                 update = 1000, missing = NA, rfit = TRUE, sign_check = FALSE, sign_eps = -.1, rs = FALSE,
                 auto_stop=FALSE,max_conv=10, rseed = 12345, digits = 4, alas = FALSE, verbose = FALSE,
                 ort.fac = 0,std=TRUE,regphi=NULL,regpsx=NULL) {
  
  Q <- as.matrix(Q)
  if (nrow(Q) != ncol(dat))
    stop("The numbers of items in data and Q are unequal.", call. = FALSE)
  
  if (iter == 0)
    stop("Parameter iter must be larger than zero.", call. = FALSE)
  
  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
    set.seed(rseed)
    
    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)
    
    Y <- t(dat)
    Y[which(Y == missing)] <- NA
    ysig<-apply(Y, 1, sd, na.rm=TRUE)
    ybar<-apply(Y, 1, mean, na.rm=TRUE)
    
    N <- ncol(Y)
    J <- nrow(Y)
    if(std){int<-FALSE}else{int<-TRUE}
    #int<-F #intercept retained or not
    
    K <- ncol(Q)
    Jp <- length(cati)
    if (Jp == 1 && cati == -1) {
      cati <- c(1:J)
      Jp <- J
    }
    
    Nmis <- sum(is.na(Y))
    mind <- which(is.na(Y), arr.ind = TRUE)
    
    const <- list(N = N, J = J, K = K, Q = Q, cati = cati, Jp = Jp, Nmis = Nmis, cand_thd = cand_thd, int = int)
    
    ######## Init ########################################################
    miter <- iter + burn
    mOmega <- array(0, dim = c(K, N))  # mean latent variable omega, K*N
    NLA <- sum(Q != 0)  #number of all lambda need to be estimated
    ELA <- array(0, dim = c(iter, NLA))  #Store retained trace of Lambda
    EMU <- array(0, dim = c(iter, J))  #Store retained trace of MU
    # EPSX <- array(0, dim = c(iter, J, J)) #Store retained trace of PSX EPHI <- array(0, dim = c(iter,
    # K, K)) #Store retained trace of PHI
    if (LD) {
      EPSX <- array(0, dim = c(iter, J * (J + 1)/2))
    } else {
      EPSX <- array(0, dim = c(iter, J))
    }
    # EPSX <- ifelse(LD,array(0, dim = c(iter, J*(J+1)/2)),array(0, dim = c(iter, J))) Store retained
    # trace of PSX
    EPHI <- array(0, dim = c(iter, K * (K - 1)/2))  #Store retained trace of PHI
    Egammas <- array(0, dim = c(iter, 1))  #Store retained trace of shrink par gammas
    Egammal <- array(0, dim = c(iter, K))  #Store retained trace of shrink par gammal (per factor)
    # Delta<<-array(0,dim=c(iter,J)) #Store retained trace of Ys scale
    Eppmc <- NULL
    if (PPMC)
      Eppmc <- array(0, dim = c(iter))
    
    init <- init(y = Y, const = const)
    Y <- init$y
    const <- init$const
    prior <- init$prior
    PSX <- init$PSX
    inv.PSX <- chol2inv(chol(PSX))
    # PHI <- init$PHI
    LA <- init$LA
    THD <- init$THD
    gammas <- init$gammas
    gammal_sq <- init$gammal_sq
    
    if (Jp > 0) {
      mnoc <- const$mnoc
      Etd <- array(0, dim = c(iter, Jp, mnoc - 1))
    }
    accrate <- 0
    if (Nmis > 0)
      Y[mind] <- rnorm(Nmis)
    
    yps<-0
    
    # OME <- t(mvrnorm(N,mu=rep(0,K),Sigma=diag(1,K))) # J*N
    sign_sw <- rep(0, K)
    # sign_eps <- -.5
    
    Eigen <- array(0, dim = c(iter, K))  #Store retained trace of Eigen
    tmp <- which(Q!=0,arr.ind=TRUE)
    pof <- matrix(0,NLA,K) #pos of est lam for each factor
    for (k in 1:K){
      ind<-(tmp[,2]==k)
      pof[ind,k]<-1
    }
    lsum <- 0
    sy <- 0
    no_conv <- 0
    LA_OF <- .99 #overflow value
    overf <- 0
    
    PHI <- diag(K)
    orl<-length(ort.fac)
    if (orl > 1){
      if (orl != K)
        stop("ort.fac should be either scalar 0/1, or a vector of 0/1 with length K.", call. = FALSE)
      nort.K <- sum(!ort.fac)
      s_PHI<-prior$s_PHI[!ort.fac,!ort.fac]
    }
    
    
    inv.PHIcov<-chol2inv(chol(PHI))
    PHIcov<-PHI
    ### for Lasso phi
    indmx <- matrix(1:K^2, nrow =K, ncol = K)
    temp <- indmx[upper.tri(indmx)]
    upperind <- temp[temp > 0]
    indmx_t <- t(indmx)
    temp <- indmx_t[upper.tri(indmx_t)]
    lowerind <- temp[temp > 0]
    const$upperind_phi <- upperind
    const$lowerind_phi <- lowerind
    ind_nod <- array(0, dim = c(K - 1, K))
    for (k in 1:K) {
      if (k == 1) {
        tmp <- 2:K
      } else if (k == K) {
        tmp <- 1:(K - 1)
      } else tmp <- c(1:(k - 1), (k + 1):K)
      ind_nod[, k] <- tmp
    }  # end of k
    const$ind_nod_phi <- ind_nod
    
    ## for horse phi
    Lambda_sq <- matrix(1, K, K)
    Nu<- matrix(1, K, K)
    xi=1
    
    ## for horse psx
    Lambda_sq_psx <- matrix(1, J, J)
    Nu_psx<- matrix(1, J, J)
    xi_psx=1
    
    ##### for ssp phi
    # Create adjacency matrix based on C
    #adj <- abs(PHI) > 1e-5
    adj <- abs(PHI) >=0 
    #adj[1:3,1:3]=T
    v0 = 0.02^2; 
    h = 50^2; 
    v1 = h * v0;
    pii <- 3/ (K- 1)  #3
    if(K<4){pii<-0.3}
    V0 <- matrix(v0, nrow = K, ncol = K)
    V1 <- matrix(v1, nrow = K, ncol = K)
    taus <- V0
    taus[adj] <- v1
    
    ##### for ssp psx
    adj_psx <- abs(PSX) > 1e-5 
    v0_psx = 0.1^2; 
    h_psx = 10^2; 
    v1_psx = h_psx * v0_psx;
    pii_psx <- 3 / (J- 1)
    V0_psx <- matrix(v0_psx, nrow = J, ncol = J)
    V1_psx <- matrix(v1_psx, nrow = J, ncol = J)
    taus_psx <- V0_psx
    taus_psx[adj_psx] <- v1_psx
    ######## end of Init #################################################
    
    ptm <- proc.time()
    for (ii in 1:miter) {
      # i <- 1
      g = ii - burn
      OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
      # PHI<-MH_PHI(ph0=PHI,ome=Omega)
      
      if (LD) {
        if(is.null(regpsx)){
          tmp <- Gibbs_PSX(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, const = const,
                           prior = prior)
          PSX <- tmp$obj
          inv.PSX <- tmp$inv
          gammas <- tmp$gammas
        }else if(regpsx=='ssp'){
          tmp <- Gibbs_PSX_ssp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, const = const,
                               pii=pii_psx,taus=taus_psx,lambda=1,V0=V0_psx,V1=V1_psx)
          PSX <- tmp$obj
          inv.PSX <- tmp$inv
          taus_psx <- tmp$taus
        }else if(regpsx=='horse'){
          tmp <- Gibbs_PSX_horse(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, const = const,
                                 Lambda_sq = Lambda_sq_psx,Nu=Nu_psx,xi=xi_psx)
          PSX <- tmp$obj
          inv.PSX <- tmp$inv
          Lambda_sq_psx = tmp$Lambda_sq
          Nu_psx=tmp$Nu
          xi_psx=tmp$xi
          }
        
        LAY <- GwMH_LA_MYC(y = Y, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                           const = const, prior = prior, alas = alas)
      } else {
        LAY <- GwMH_LA_MYE(y = Y, ome = OME, la = LA, psx = PSX, gammal_sq = gammal_sq, thd = THD,
                           const = const, prior = prior, alas = alas)
        PSX <- LAY$psx
        inv.PSX <- chol2inv(chol(PSX))
      }
      
      LA <- LAY$la # OME <- LAY$ome
      if(std){
        tmp1<-abs(LA)>LA_OF
        if (any(tmp1)){
          overf <- overf+(colSums(tmp1)>0)
          # LA1[tmp1]<-LA[tmp1]
          LA[LA>LA_OF]<-LA_OF
          LA[LA< -LA_OF]<--LA_OF
        }
      }
      
      
      if (sign_check) {
        chg <- (colSums(LA)<= sign_eps)
        # chg <- (colSums(LA)<= LA_eps)
        if (any(chg)) {
          sign <- diag(1 - 2 * chg)
          sign_sw <- sign_sw + chg
          # if(g<0){chg0_count <- chg0_count + chg}else{chg_count <- chg_count + chg}
          LA <- LA %*% sign
          OME <- t(t(OME) %*% sign)
          if(verbose){
            print(c("ii=", ii), quote = FALSE)
            cat(sign_sw, fill = TRUE, labels = "#Sign switch:")
          }
        }
      } #end if
      
      gammal_sq <- LAY$gammal_sq
      # OME <- Gibbs_Omega(y = Y, la = LA, phi = PHI, inv.psx = inv.PSX, N = N, K = K)
      
      if(is.null(regphi)){
        ### original phi
        if (orl > 1){
          # phi0<-PHI[!ort.fac,!ort.fac]
          tmp <- MH_PHI(phi = PHI[!ort.fac,!ort.fac], ome = OME[!ort.fac,], N = N, K = nort.K, s0 = s_PHI,std=std)
          PHI[!ort.fac,!ort.fac]<-tmp
        }else if (ort.fac == 0){
          PHI <- MH_PHI(phi = PHI, ome = OME, N = N, K = K, s0 = prior$s_PHI,std=std)
        }
      }else if(regphi=='lasso'){
        ### lasso phi
        if (orl > 1){
          # phi0<-PHI[!ort.fac,!ort.fac]
          tmp <- MH_PHI_LASSO(phi = PHIcov[!ort.fac,!ort.fac], ome = OME[!ort.fac,],inv.phi= inv.PHIcov[!ort.fac,!ort.fac],const=const, prior=prior)
          PHIcov[!ort.fac,!ort.fac]<-tmp$phi
          inv.PHIcov[!ort.fac,!ort.fac] <- tmp$inv.phi
        }else if (ort.fac == 0){
          tmp <- MH_PHI_LASSO(phi = PHIcov, ome = OME,inv.phi= inv.PHIcov,const=const, prior=prior)
          PHIcov<-tmp$phi
          inv.PHIcov<- tmp$inv.phi
        }
        if(std){
          tmp <- chol2inv(chol(sqrt(diag(diag(PHIcov)))))
          PHI <- tmp %*% PHIcov %*% tmp
        }else{PHI<-PHIcov}
        
      }else if(regphi=='horse'){
        #### horse phi
        if (orl > 1){
          # phi0<-PHI[!ort.fac,!ort.fac]
          tmp <- MH_PHI_horse(phi = PHIcov[!ort.fac,!ort.fac], ome = OME[!ort.fac,],inv.phi= inv.PHIcov[!ort.fac,!ort.fac],const=const, Lambda_sq = Lambda_sq[!ort.fac,!ort.fac],Nu=Nu[!ort.fac,!ort.fac],xi=xi)
          PHIcov[!ort.fac,!ort.fac]<-tmp$phi
          inv.PHIcov[!ort.fac,!ort.fac] <- tmp$inv.phi
          Lambda_sq[!ort.fac,!ort.fac] = tmp$Lambda_sq
          Nu[!ort.fac,!ort.fac]=tmp$Nu
          xi=tmp$xi
        }else if (ort.fac == 0){
          tmp <- MH_PHI_horse(phi = PHIcov, ome = OME,inv.phi= inv.PHIcov,const=const,  Lambda_sq = Lambda_sq,Nu=Nu,xi=xi)
          PHIcov<-tmp$phi
          inv.PHIcov<- tmp$inv.phi
          Lambda_sq = tmp$Lambda_sq
          Nu=tmp$Nu
          xi=tmp$xi
        }
        if(std){
          tmp <- chol2inv(chol(sqrt(diag(diag(PHIcov)))))
          PHI <- tmp %*% PHIcov %*% tmp
        }else{PHI<-PHIcov}
        
      }else if(regphi=='ssp'){
        #### horse phi
        if (orl > 1){
          # phi0<-PHI[!ort.fac,!ort.fac]
          tmp <- MH_PHI_ssp(phi = PHIcov[!ort.fac,!ort.fac], ome = OME[!ort.fac,],inv.phi= inv.PHIcov[!ort.fac,!ort.fac],const=const, pii=pii,taus=taus[!ort.fac,!ort.fac],lambda=1,V0=V0[!ort.fac,!ort.fac],V1=V1[!ort.fac,!ort.fac])
          PHIcov[!ort.fac,!ort.fac]<-tmp$phi
          inv.PHIcov[!ort.fac,!ort.fac] <- tmp$inv.phi
          taus[!ort.fac,!ort.fac]<-tmp$taus
        }else if (ort.fac == 0){
          tmp <- MH_PHI_ssp(phi = PHIcov, ome = OME,inv.phi= inv.PHIcov,const=const,pii=pii,taus=taus,lambda=1,V0=V0,V1=V1)
          PHIcov<-tmp$phi
          inv.PHIcov <- tmp$inv.phi
          taus<-tmp$taus
        }
        if(std){
          tmp <- chol2inv(chol(sqrt(diag(diag(PHIcov)))))
          PHI <- tmp %*% PHIcov %*% tmp
        }else{PHI<-PHIcov}
        
      }

      
      if (Jp > 0) {
        Y[cati, ] <- LAY$ys
        THD <- LAY$thd
        accrate <- accrate + LAY$accr
      }
      if (Nmis > 0)
        Y[mind] <- LAY$ysm[mind]
      
      # Save results
      if ((g > 0)) {
        mOmega <- mOmega + OME
        ELA[g, ] <- LA[Q != 0]
        Eigen[g,] <- (ELA[g,]^2)%*%pof
        # EPSX[g, , ] <- PSX EPSX[g, ] <- ifelse(LD,PSX[lower.tri(PSX,diag=TRUE)],diag(PSX))
        if (LD) {
          EPSX[g, ] <- PSX[lower.tri(PSX, diag = TRUE)]
        } else {
          EPSX[g, ] <- diag(PSX)
        }
        Egammas[g, ] <- gammas
        Egammal[g, ] <- colMeans(sqrt(gammal_sq))
        # EPHI[g, , ] <- PHI[, ]
        EPHI[g, ] <- PHI[lower.tri(PHI)]
        # EMU[g,]<-MU
        if (Jp > 0)
          Etd[g, , ] <- THD[, 2:mnoc]
        if (PPMC)
          Eppmc[g] <- post_pp(y = Y, ome = OME, la = LA, psx = PSX, inv.psx = inv.PSX, N = N,
                              J = J)
        
        if (rfit){
          Yc <- Y - LA %*% OME  # J*N
          tmp<-(t(Yc) %*%chol(inv.PSX))^2
          # tmp<-(t(Yc) %*%chol(chol2inv(chol(PSX))))^2
          lsum<-lsum+sum(tmp)+N*(log(det(PSX)))
        } #end dic
        
        if (rs) {
          # if (Nmis == 0){
          #     ytmp<-LAY$ysm
          #   }else{
          #     ytmp<-Y[mind]
          #   }
          yps<-yps + LAY$ysm
        }
      } #end g
      
      if (ii%%update == 0){
        if (g > 0) {
          
          APSR <- schain.grd(Eigen[1:g,])
          # if (auto_stop) {
          if (max(APSR[,1]) < 1.1) {
            no_conv <- no_conv + 1
          } else{
            no_conv <- 0
          }
          # } # end auto_stop
        } # end g
        
        
        if(verbose){
          # print(proc.time() - ptm)
          Shrink <- colMeans(sqrt(gammal_sq))
          Feigen <- diag(crossprod(LA))
          NLA_le3 <- colSums(abs(LA) >= 0.3)
          # Meigen <- colMeans(Eigen)
          # Mlambda<-colMeans(LA)
          
          cat(ii, fill = TRUE, labels = "\nTot. Iter =")
          # print(rbind(Feigen, NLA_le3, Shrink,sign_sw))
          print(rbind(Feigen, NLA_le3, Shrink))
          cat(overf, fill = TRUE, labels = 'Overflow')
          if (g > 0) cat(c(t(APSR[,1]),no_conv), fill = TRUE, labels = "EPSR & NCONV")
          
          if (Jp > 0) {
            cat(colMeans(THD), fill = TRUE, labels = "Ave. Thd:")
            cat(accrate[1]/update, fill = TRUE, labels = "Acc Rate:")
            accrate <- 0
          }
          
          if (LD) {
            opsx <- abs(PSX[row(PSX) != col(PSX)])
            tmp <- c(sum(opsx > 0.2)/2, sum(opsx > 0.1)/2, gammas)
            print(c("LD>.2", ">.1", "Shrink"), quote = FALSE)
            print(tmp)
          }
          
        }#end verbose
        
        if (auto_stop * no_conv >= max_conv) break
      }  # end update
      
    }  #end of g MCMAX
    
    if(verbose){
      # cat(chg0_count,chg_count, fill = TRUE, labels = "\n#Sign change:")
      print(proc.time()-ptm)
    }
    
    # if (conv != 0) {
    #     st <- g/2 + 1
    #     ELA <- ELA[(st):g, ]
    #     # EPSX <- EPSX[(st):g, , ] EPHI <- EPHI[(st):g, , ]
    #     EPSX <- EPSX[(st):g, ]
    #     EPHI <- EPHI[(st):g, ]
    #     Egammal <- Egammal[(st):g, ]
    #     Egammas <- Egammas[(st):g, ]
    #     if (Jp > 0)
    #         Etd <- Etd[(st):g, , ]
    #     Eppmc <- Eppmc[(st):g]
    #     mOmega <- mOmega/2
    #     iter <- burn <- g/2
    # }
    
    if (auto_stop * no_conv >= max_conv) {
      ELA <- ELA[1:g, ]
      Eigen <- Eigen[1:g, ]
      EPSX <- EPSX[1:g, ]
      EPHI <- EPHI[1:g, ]
      Egammal <- Egammal[1:g, ]
      Eppmc <- Eppmc[1:g]
      Egammas <- Egammas[1:g, ]
      if (Jp > 0)
        Etd <- Etd[1:g, , ]
      iter <- g
    }
    
    if (rfit){
      lpry <- lsum/iter+N*log(2*pi)
      # lsum <- lsum/iter
    }
    
    # chg1_count<-rbind(chg0_count,chg_count)
    out <- list(Q = Q, LD = LD, LA = ELA, Omega = mOmega/iter, PSX = EPSX, iter = iter, burn = burn,
                PHI = EPHI, gammal = Egammal, gammas = Egammas, Nmis = Nmis, PPP = Eppmc, Eigen = Eigen,
                auto_conv = c(auto_stop, no_conv, max_conv), Y = Y, lpry=lpry, time = (proc.time()-ptm))
    
    if (Jp > 0) {
      out$cati = cati
      out$THD = Etd
      out$mnoc = mnoc
    }
    
    if (rs) {
      yp<-yps/iter*ysig+ybar
      out$yp = t(yp)
    }
    
    class(out) <- c("lawbl")
    
    if (!is.null(oldseed))
      .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)
    
    # options(digits = old_digits)
    
    return(out)
}

