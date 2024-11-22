#### test reg phi
rm(list=ls())
setwd("C:/Users/zhangyf/Downloads/LAWBL/R")
source("sim_lvm.R")
source("utils.R")
source("Init.R")
source("Omega_Phi.R")
source("Gibbs_PSX.R")
source("GwMH_LA_MYE.R", local = TRUE)
source("GwMH_LA_MYC.R", local = TRUE)
source("pcfa_regphi.R")

library(MASS)
library(coda)
library(stats)
library(LAWBL)


################ Non-standardized
library(lavaan)
population.model <- ' f1 =~ x1 + 0.7*x2 + 0.7*x3 + 0.7*x4 + 0.7*x5 + 0.7*x6 
                      f2 =~ x7 + 0.7*x8 + 0.7*x9 + 0.7*x10 + 0.7*x11 + 0.7*x12 
                      f3 =~ x13 + 0.7*x14 + 0.7*x15 + 0.7*x16 + 0.7*x17 + 0.7*x18 
                      f1 ~ 0.3*f2 
                      f3 ~ 0.3*f1 + 0.3*f2
                    '
set.seed(1234)
y <- simulateData(population.model, sample.nobs=1000L)
K=3;J=18;
Q<-matrix(0,J,K)
Q[1:6,1]<-Q[7:12,2]<-Q[13:18,3]<-1
Q[1,1]<-Q[7,2]<-Q[13,3]<-999  ##fix loading to 1
m0 <- pcfa(dat = y, Q = Q, LD =F,burn =2000, iter = 2000,sign_check = T,verbose = T,alas=T)
m1 <- pcfa_regphi(dat = y,std=F, Q = Q,regphi = 'lasso',regpsx = NULL, LD = F,burn = 2000, iter = 2000,sign_check = F,verbose = T,alas=T)
summary(m0, what = 'qlambda',detail=T)
summary(m0, what = 'phi',detail=F)

################ Standardized
y <- sim18cfa0$dat
Q<-matrix(-1,J,K)
Q[1:6,1]<-Q[7:12,2]<-Q[13:18,3]<-1
m2 <- pcfa_regphi(dat = y,std=T, Q = Q,regphi = 'ssp',regpsx = NULL, LD = F,burn = 2000, iter = 2000,sign_check = F,verbose = T,alas=T)
summary(m2, what = 'qlambda',detail=T)
summary(m2, what = 'phi',detail=F)


K=3;J=18;N=1000;lac=0;cpf=0
mla=0.7;PHI<-0.3; LD<-F;estLD=F;cati=NULL;noc=NULL
if (LD==T){ecr<-.2;necw<-2;necb<-2} else {ecr<-0;necw<-0;necb<-0}
Q<-matrix(-1,J,K)
out <- sim_lvm(N=N,lam=mla,J=J,K=K,lac = lac,cpf=cpf,phi=PHI,rseed=123,ecr=ecr,necw=necw,necb=necb)
y<-out$dat
m3<- pcfa_regphi(dat = y,std=T, Q = Q,regphi = 'ssp',regpsx = NULL, LD = F,burn = 2000, iter = 2000,sign_check = F,verbose = T,alas=T)
