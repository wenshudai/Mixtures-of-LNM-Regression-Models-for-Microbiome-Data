rm(list=ls())

library(curatedMetagenomicData)
library(microbiome)
library(ALDEx2)
library(mia)
library(mclust)
library(MASS)
library(nnet)
library("plyr")
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# setwd("G:/My Drive/LNMVGA/")
source("/Users/alexia/Documents/Bing/Project2/code/Real data/hmp2012/t_functions.R")
rawdata <- read.delim("/Users/alexia/Documents/Bing/Project2/code/Real data/jrssc/real data/Phylum.txt",header = T, sep = " ")

#################################################################
data <- rawdata[c(2,3,4)][1:259,]
data[data$Phenotype==1,1] <-5
data[data$Phenotype==4,1] <-1
data[data$Phenotype==5,1] <-4

rawdata <- rawdata[1:259,]
output <- data.matrix(rawdata)

most_ab<- HMP::Data.filter(output[,-(1:5)][1:259,], "data", 0, 5)
colnames(most_ab)

# X <- cbind(1,factor(data$Phenotype),factor(data$Gender),data$Standardized.Age)
X.datframe <- data.frame(int=1,Pheno=factor(data$Phenotype),Gender=factor(data$Gender),Age=data$Standardized.Age)
## factor 
X.f <- cbind(1,model.matrix(~ X.datframe$Pheno)[,-1],X.datframe$Gender,data$Standardized.Age)
## numeric 
X.n <- cbind(1,X.datframe$Pheno,model.matrix(~ X.datframe$Gender)[,-1],Age=data$Standardized.Age)
W <- most_ab
P=ncol(X.n)

# W_init <- W

n = nrow(W)
K=ncol(W)-1

# W <- W_count
W_init <- W
if(any(W_init==0)){W_init[which(W_init==0)]<-1}

n = nrow(W)
K=ncol(W)-1

multi_const<<-lgamma(rowSums(W)+1)-apply(lgamma(W+1),1,sum)
W_init <- W
if(any(W_init==0)) W_init[which(W_init==0)]<-1
n <- nrow(W)
K<-ncol(W)-1 

# true_lab<-Data$label

theta <- W_init/rowSums(W_init)
latent <- log(theta/theta[,(K+1)])

all_G <- c(1:5)

loglik_list <- vector("list",5)
par_stor_list <- vector("list",5)

BIC <- numeric(5)
ICL <- numeric(5)
ARI <- numeric(5)

for(G in all_G){
  # G = 2
  ###Initialization
  z <- unmap(kmeans(latent,center=G,nstart=50)$cluster)
  pi_g <- colSums(z)/sum(z)
  xi_hat <- rep(1,n)
  m_hat <- latent
  V_hat <- vector("list",n)
  Vmat<-matrix(rep(c(rep(0.1,K),0),n),nrow=n,byrow=T)
  
  for(i in 1:n){V_hat[[i]]<-diag(Vmat[i,])}
  
  nu_hat <- rep(50,G)
  u_hat <- matrix(1,n,G)
  
  musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=z,u=u_hat)
  if(all(!is.na(musig_up$mu))){
    mu_hat <- musig_up$mu
    Sig_hat <- musig_up$Sig
  }else{
    mu_hat <- matrix(rep(colMeans(m_hat),G),nrow=G,ncol=(K+1),byrow=TRUE)
    Sig_hat <- list()
    for(g in 1:G){Sig_hat[[g]] <- rbind(cbind(diag(1,K),0),0)}
  }
  red_sig<-list()
  for(g in 1:G){
    red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
  }
  iSig_hat<-list()
  for(g in 1:G){
    iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
  }
  
  s_hat <- (K+1+nu_hat)/2
  r_hat <- matrix(0,n,G)
  for (g in 1:G) {
    r_hat[,g] <- (mahalanobis(m_hat, mu_hat[g,],cov=iSig_hat[[g]],inverted=TRUE)+nu_hat[g])/2 
  }
  
  # zhat <- z
  
  #######################EM loop starts######################
  ###########################################################
  aloglik<-NULL
  loglik<-NULL
  aloglik[c(1,2,3,4,5)]<-0
  it_max<-500
  
  
  it <- 2
  stop <- 0
  while(stop < 1){
    condition <- 0
    second_vec <- numeric(G)
    # first_vec <- numeric(G)
    for_hic <- matrix(NA,nrow=n,ncol=G)
    ##Updating ELBO
    temp_ELBO<-matrix(NA,nrow=n,ncol=G)
    for (g in 1:G){
      first_e<-diag(W%*%t(m_hat))
      second_e<-log(rowSums(exp(m_hat+Vmat/2)))*rowSums(W)
      for_hic[,g] <- first_e-second_e+multi_const
      # second_e<-(rowSums(exp(m_hat+Vmat/2))/xi_hat+log(xi_hat)-1)*rowSums(W)
      third_e<-0.5*log(det(red_sig[[g]]))
      fourth_e<-0.5*mahalanobis(m_hat,center=mu_hat[g,],cov=iSig_hat[[g]],inverted=TRUE)
      fifth_fun<-function(V_hatm){
        sum(diag(iSig_hat[[g]]%*%V_hatm))
      }
      fifth_e<-0.5*sapply(V_hat,fifth_fun)
      sixth_e<-0.5*rowSums(log(Vmat[,-(K+1)]))
      srnu_e<-s_hat[g]-log(r_hat[,g])+log(gamma(s_hat[g]))-(s_hat[g]-1)*digamma(s_hat[g])+(K+nu_hat[g]-1)*0.5*(digamma(s_hat[g])-log(r_hat[,g]))-0.5*nu_hat[g]*u_hat[,g]+0.5*nu_hat[g]*log(0.5*nu_hat[g])-log(gamma(0.5*nu_hat[g]))
      temp_ELBO[,g]<-first_e-second_e-third_e-fourth_e*u_hat[,g]-fifth_e*(K+1)*u_hat[,g]+sixth_e+(K+1)/2+multi_const+srnu_e
      
    }
    pi_mat<-matrix(pi_g,ncol=G,nrow=n,byrow=TRUE)
    for_tau<-pi_mat*exp(temp_ELBO)
    for_ll<-rowSums(for_tau)
    if (it<5){
      for_ll[for_ll==0]<-1
      loglik[it]<-sum(log(for_ll))
      for_z<-pi_mat*exp(temp_ELBO-apply(temp_ELBO,1,max))
      tau<-for_z/rowSums(for_z)} else {
        if(any(for_ll==0)){
          stop <- 1
          condition <- 1
          loglik[it] <- NA
          it<-it+1
        }
        loglik[it]<-sum(log(for_ll))
        tau<-for_tau/rowSums(for_tau)
      }
    zhat<-tau
    
    if(condition < 1){
      ########Update variational parameters
      varpar_up <- try(varpar_fun(W=W,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,K=K,z=zhat,iSig=iSig_hat,u=u_hat),silent=TRUE)
      if(class(varpar_up)=="try-error"){
        xi_hat <- NA
        m_hat <- NA
        V_hat <- NA
        Vmat<-NA
        varpar_up <- list(xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat)
      }else{
        xi_hat <- varpar_up$xi
        m_hat <- varpar_up$m
        V_hat <- varpar_up$V
        Vmat<-varpar_up$Vmat
      }
      
      if(any(is.na(varpar_up))){
        stop<-1
        loglik[it] <- NA
        it<-it+1
      }else{
        ####    # M-step
        pi_g <- colSums(zhat)/sum(zhat)
        ng<-colSums(zhat)
        musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=zhat,u=u_hat)
        mu_hat <- musig_up$mu
        Sig_hat <- musig_up$Sig
        # mu_hat <- cbind(rbind(true_mu[[1]],true_mu[[2]]),0)
        # Sig_hat <- true_sig
        # Sig_hat[[1]] <- rbind(cbind(Sig_hat[[1]],0),0)
        # Sig_hat[[2]] <- rbind(cbind(Sig_hat[[2]],0),0)
        
        iSig_hat<-list()
        for(g in 1:G){
          iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
        }
        red_sig<-list()
        for(g in 1:G){
          red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
        }
        
        ######## update s,r,nu,u
        nu_hat<-rep(0,G)
        for (g in 1:G){
          temp<-function(v){digamma(s_hat[g])-1/ng[g]*sum(zhat[,g]*((log(r_hat[,g])+s_hat[g]/r_hat[,g])))+1-digamma(v/2)+log(v/2)}
          nu_hat[g]<-try(uniroot(temp, lower=0.0001, upper=1000)$root,silent=TRUE)
          if (nu_hat[g]>200|class(nu_hat[g])=="try-error") nu_hat[g]<-200
        }
        s_hat <- (K+1+nu_hat)/2
        r_hat <- matrix(0,n,G)
        for (g in 1:G) {
          r_hat[,g] <- (mahalanobis(m_hat, mu_hat[g,],cov=iSig_hat[[g]],inverted=TRUE)+nu_hat[g])/2 
        }
        if(G == 1){u_hat <- s_hat/r_hat}else{u_hat <- t(apply(r_hat,1,function(x) s_hat/x))}
        
        # nu_hat <- rep(200,G)
        
        if (it>5){
          #Aitkaine's stopping criterion
          if ((loglik[it-1]-loglik[it-2])==0) checks<-1 else{
            a<-(loglik[it]-loglik[it-1])/(loglik[it-1]-loglik[it-2])
            add_to<-(1/(1-a)*(loglik[it]-loglik[it-1]))
            # }
            aloglik[it]<-loglik[it-1]+add_to
            if (abs(aloglik[it]-aloglik[it-1])<0.001) stop<-1 else stop<-stop
          }
        }	
        cat("G =",G,", iteration =",it,"\n")
        it<-it+1
        if (it==it_max) stop<-1
      }
    }
  }
  
  npar <- (K+1)*G+0.5*(K+1)*(K)*G+(G-1)
  BIC[G]<-2*loglik[(it-1)]-npar*log(n)
  mapz<-mclust::unmap(mclust::map(zhat))
  forICL<-function(g){sum(log(zhat[which(mapz[,g]==1),g]))}
  ICL[G]<-BIC[G]+2*sum(sapply(1:ncol(mapz),forICL))
  # ARI[G]<-adjustedRandIndex(mclust::map(zhat),true_lab)
  loglik_list[[G]] <- loglik
  par_stor_list[[G]] <- list(m=m_hat,V=V_hat,Vm=Vmat,xi=xi_hat,mu=mu_hat,Sig=Sig_hat,iSig=iSig_hat,pi_g=pi_g,nu=nu_hat,s=s_hat,r=r_hat,u=u_hat,z=zhat,class=mclust::map(zhat))
}
which.max(BIC)

save.image("/Users/alexia/Documents/Bing/Project2/code/Real data/jrssc/j_214.Rdata")



