library(expm)
library(MASS)
construct_B<-function(x,knots)
{
  #Correlation
  Omega_Z0=matrix(unlist(lapply(
    knots,FUN=function(x){abs((x-knots)^3)}
  )),nrow=length(knots),byrow=T)
  diag(Omega_Z0)=1
  return(cbind(1, x, matrix(unlist(lapply(
    x,
    FUN = function(x) {
      abs((x - knots) ^ 3)
    }
    #  )), nrow = length(x), byrow = T)%*%sqrtm(solve(Omega_Z0)) ))
  )), nrow = length(x), byrow = T) ))
}

#time_grids_obs=unique(Complete.data$t)

knots=quantile(Complete.data$t,seq(0,1,length=l+2))[-c(1,l+2)]

B_list=aggregate(Complete.data$t,by=list(Complete.data$sname),FUN=
                   function(x){construct_B(x,knots)})
#B_list=aggregate(Complete.data$t,by=list(Complete.data$sname),FUN=
#                   function(x){construct_B_poly(x,1:l)})
MB_matrix=B_list$x
Time_Grid=seq(0,1,length=Grid.num)

#Construct Time Grid B
MB_Grid=construct_B(Time_Grid,knots)
#MB_Grid=construct_B_poly(Time_Grid,1:l)
MB_Grid_Product=t(MB_Grid)%*%MB_Grid


n_unit=length(unique(Complete.data$sname))
n_total=nrow(Complete.data)

##Cluster ID
# cluster_label_value=unique(Complete.data$grp)
# cluster_label=1:length(cluster_label_value)
index_cluster=aggregate(1:nrow(Complete.data),by=list(Complete.data$grp),FUN=function(x){x})
#12 cluster in total
##index label for each cluster
cluster_label=unlist(lapply(Complete.data$grp,FUN=function(x){which(index_cluster$Group.1==x)}))
##The index set for each cluster
index_cluster=index_cluster$x

#Hydro
hydro_cluster=aggregate(1:nrow(Complete.data),by=list(Complete.data$hydroyear),FUN=function(x){x})
hydro_cluster_label=unlist(lapply(Complete.data$hydroyear,FUN=function(x){which(hydro_cluster$Group.1==x)}))
##The index set for each cluster
hydro_cluster=hydro_cluster$x

#Covariate Matrix
#covariate.dsi.index="n_adults"
#X_M=as.matrix(Complete.data[,covariate.index.dsi])
X_M = model.matrix(as.formula(paste(
  "~", paste(covariate.index.dsi, collapse = "+")
)), Complete.data)
X_M=X_M[,-1]

#Centralize
X_M=scale(X_M)
X_M_mean=attr(X_M,"scaled:center")
X_M_sd=attr(X_M,"scaled:scale")
if(simple==1){
  X_M=cbind(rep(1,nrow(X_M)),X_M)
}else{
  X_M=cbind(construct_B(Complete.data$t,knots),X_M)}
#X_M=cbind(construct_B(Complete.data$t,knots),X_M,Complete.data$DSI_F)
X_Product=t(X_M)%*%X_M


id_index=aggregate(1:nrow(X_M),by=list(Complete.data$sname),FUN=function(x){x})
id_index=id_index$x




#Outcome
Obs=Complete.data[,mediation.index[m]]
Obs_list=aggregate(Complete.data[,mediation.index[m]],by=
                     list(Complete.data$sname),
                   function(x){x})
# Obs=Complete.data$gc
# Obs_list=aggregate(Complete.data$gc,by=
#                      list(Complete.data$sname),
#                    function(x){x})
Obs_value=Obs_list$x

#Treatment Indicator
Treatment_Ind=aggregate(Complete.data[,adverse.index[adv_id]],by=list(Complete.data$sname),
                        FUN=function(x){if(any(x)){return(1)}else{return(0)}})$x
N_0=sum(Treatment_Ind==0)
N_1=n_unit-N_0

#Penalty Function
Omega=diag(c(0,0,rep(1,l)))

#MCMC Section, Preparation
#Parameter for Covariates
theta=matrix(0,ncol=ncol(X_M),nrow=simu_time)
#Principal Score
Beta_ini=matrix(runif(K*n_unit,min=seq(-2,-0.2,length=K),
                      max=seq(2,0.2,length=K)),ncol=K,nrow=n_unit,
                byrow=T)
Beta=rep(list(NA),simu_time)
Beta[[1]]=Beta_ini

#Precision Parameter
lambda=matrix(2,ncol=K,nrow=simu_time)
Psi_ini=matrix(runif(K*(l+2),-1,1),ncol=K,nrow=l+2)
Psi=rep(list(NA),simu_time)
Psi[[1]]=Psi_ini


#Precison Parameters
sigma_noise=rep(0,simu_time);sigma_noise[1]=1
sigma_K=matrix(0,ncol=K,nrow=simu_time);sigma_K[1,]=seq(2,1,length=K)

sigma_theta=100^2

#TE Parameter
tau_0=tau_1=matrix(0,ncol=K,nrow=simu_time)

# #Process value
# m_process_value=matrix(0,nrow=simu_time,ncol=nrow(Complete.data))

#Parameter for MGP
#For causal parameters
tau_delta=matrix(0,ncol=K,nrow=simu_time);tau_delta[1,]=rep(1,K)
#For RD
rd_delta=matrix(0,ncol=K,nrow=simu_time);rd_delta[1,]=rep(1.1,K)
tau_delta=matrix(0,ncol=K,nrow=simu_time);tau_delta[1,]=rep(1.1,K)

#rd_xi=matrix(0,ncol=n_unit,nrow=simu_time)
sigma_noise=rep(0,simu_time);sigma_noise[1]=1
kappa=matrix(1,ncol=n_unit,nrow=simu_time);

sigma_K=matrix(0,ncol=K,nrow=simu_time);
sigma_K[1,1]=1/rd_delta[1,1]
for (k in 2:K)
{
  sigma_K[1,k]=sigma_K[1,k-1]/rd_delta[1,k]
}
if(tau_mgp==1){
  sigma_tau=matrix(0,ncol=K,nrow=simu_time)
  sigma_tau[1,1]=1/tau_delta[1,1]
  for (k in 2:K)
  {
    sigma_tau[1,k]=sigma_tau[1,k-1]/tau_delta[1,k]
  }}else{
    sigma_tau=rep(1,simu_time)
  }

#Hyperparameter for Gamma
amu=ard=matrix(0,ncol=2,nrow=simu_time);amu[1,]=c(1,2)
ard=matrix(0,ncol=2,nrow=simu_time);ard[1,]=c(1,2)

phi_ini=matrix(1,ncol=K,nrow=n_unit)
phi=rep(list(NA),simu_time)
phi[[1]]=phi_ini


##Cluster Random Effect
cluster_rd=matrix(0,ncol=length(index_cluster),nrow=simu_time)
cluster_rd_sigma=rep(1,simu_time)

##Hydro Cluster Random Effect
hydro_cluster_rd=matrix(0,ncol=length(hydro_cluster),nrow=simu_time)
hydro_cluster_rd_sigma=rep(1,simu_time)

##Decompose for the Half Cauchy distribution
# cluster_rd_xi=rep(1,simu_time)
# cluster_rd_eta=rep(1,simu_time)
#Step Size for MH
step_size=2
v=10
#HP
a_rd_1=20
a_rd_2=50
#Parameter for Half Cauchy A=25


##Gibbs Sampler
for (t in 2:simu_time){
  # if(t<=simu_time/5)
  #   truncate=1
  # else
  #   truncate=0
  #Fix k
  ####################################
  #Sample Eigen function
  ####################################
  Psi[[t]]=Psi[[t-1]]
  Beta[[t]]=Beta[[t-1]]
  for (k in sample(1:K,K)){
    Q_k=Reduce('+',lapply(
      1:n_unit,
      FUN = function(x) {
        t(MB_matrix[[x]]) %*% MB_matrix[[x]] * Beta[[t]][x, k]^2*
          kappa[t-1,x]}))/sigma_noise[t-1]+lambda[t-1,k]*Omega
    
    
    l_k=Reduce('+',lapply(1:n_unit,FUN=function(i){t(MB_matrix[[i]])%*%(Obs_value[[i]]-X_M[id_index[[i]],]%*%theta[t-1,]-
                                                                          cluster_rd[t-1,cluster_label[id_index[[i]]]]-
                                                                          hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
                                                                          MB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*Beta[[t]][i,k]*kappa[t-1,i]}))/sigma_noise[t-1]
    #Add a constraint
    #C_k=t(cbind(c(1,rep(0,l+1)),PsPsi[[t]][,-k]))%*%MB_Grid_Product
    C_k=t(Psi[[t]][,-k])%*%MB_Grid_Product
    
    Q_L=t(chol(Q_k))
    l_tilde=solve(Q_L)%*%l_k
    Psi_0=solve(t(Q_L))%*%(l_tilde+rnorm(l+2))
    
    C_k=t(Psi[[t]][,-k])%*%MB_Grid_Product
    
    C_tilde=solve(Q_L)%*%t(C_k) 
    C_tilde=solve(t(Q_L))%*%C_tilde
    
    Psi[[t]][, k] = Psi_0 - C_tilde%*% 
      solve(C_k %*% C_tilde)%*%C_k%*%Psi_0
    
    
    #Norm of f
    temp_norm=as.numeric(t(Psi[[t]][, k])%*%MB_Grid_Product%*%Psi[[t]][, k])/Grid.num
    Psi[[t]][,k]=Psi[[t]][,k]/sqrt(abs(temp_norm))
    Beta[[t]][,k]=Beta[[t]][,k]*sqrt(abs(temp_norm))
    
    ##Sample Lambda
    #lambda[t,k]=rgamma(1,shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    if(shrink==1)
    {
      low_u=pgamma(100*sigma_K[t-1,k],shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    }else{
      low_u=qgamma(10^(-8),shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    }
    u=min(runif(1,low_u,1),0.99)
    lambda[t,k]=qgamma(u,shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    
  }
  ####################################
  ##Sample Principal Score: Beta
  ####################################
  for (i in sample(1:n_unit,n_unit)){
    temp_f=MB_matrix[[i]]%*%Psi[[t]]
    f_sq_sum=apply(temp_f^2,2,sum)
    temp_sigma=1/(f_sq_sum*kappa[t-1,i]/sigma_noise[t-1]+phi[[t-1]][i,]/sigma_K[t-1,])
    for (k in sample(1:K,K))
    {
      temp_mean=
        sum((Obs_value[[i]]-X_M[id_index[[i]],]%*%theta[t-1,]-cluster_rd[t-1,cluster_label[id_index[[i]]]]-
               hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
               MB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*temp_f[,k])*kappa[t-1,i]/sigma_noise[t-1]+
        (tau_0[t-1,k]*(1-Treatment_Ind[i])+tau_1[t-1,k]*Treatment_Ind[i])*phi[[t-1]][i,k]/sigma_K[t-1,k]
      Beta[[t]][i,k]=rnorm(1,mean=temp_mean*temp_sigma[k],sd=sqrt(temp_sigma[k]))
    }
  }
  
  ####################################
  #Sample Tau_0,Tau_1
  ####################################
  if(tau_mgp==1){
    tau_0_sigma=1/(apply(phi[[t-1]][Treatment_Ind==0,],2,sum)/sigma_K[t-1,]+1/sigma_tau[t-1,])
  }else{
    tau_0_sigma=1/(apply(phi[[t-1]][Treatment_Ind==0,],2,sum)/sigma_K[t-1,]+1/sigma_tau[t-1])
  }
  tau_0_mean=apply(Beta[[t]][Treatment_Ind==0,]*phi[[t-1]][Treatment_Ind==0,],2,sum)/sigma_K[t-1,]
  tau_0[t,]=rnorm(n=K,mean=tau_0_mean*tau_0_sigma,sd=sqrt(tau_0_sigma))
  
  if(tau_mgp==1){
    tau_1_sigma=1/(apply(phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]+1/sigma_tau[t-1,])
  }else{
    tau_1_sigma=1/(apply(phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]+1/sigma_tau[t-1])
  }
  #  tau_1_sigma=1/(apply(phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]+1/sigma_tau[t-1])
  tau_1_mean=apply(Beta[[t]][Treatment_Ind==1,]*phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]
  tau_1[t,]=rnorm(n=K,mean=tau_1_mean*tau_1_sigma,sd=sqrt(tau_1_sigma))
  
  ####################################
  #Sample theta
  ####################################
  X_Product=Reduce('+',lapply(1:n_unit,FUN=function(i){
    t(X_M[id_index[[i]],])%*%X_M[id_index[[i]],]*kappa[t-1,i]}))
  
  Theta_Sigma=solve(X_Product/sigma_noise[t-1]+diag(rep(1/sigma_theta,ncol(X_M))))
  
  process_value=unlist(lapply(1:n_unit,FUN=function(x){MB_matrix[[x]]%*%Psi[[t]]%*%Beta[[t]][x,]}))
  #  Theta_Mean=Theta_Sigma%*%(t(X_M)%*%(Obs-process_value)%*%diag(kappa[t-1,]/sigma_noise[t-1]))
  Theta_Mean=Theta_Sigma%*%(
    Reduce('+',lapply(1:n_unit,
                      FUN=function(i){t(X_M[id_index[[i]],])%*%
                          (Obs[id_index[[i]]]-process_value[id_index[[i]]]-cluster_rd[t-1,cluster_label[id_index[[i]]]]-
                             hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]])*kappa[t-1,i]}))/sigma_noise[t-1])
  theta[t,]=mvrnorm(n=1,mu=Theta_Mean,Sigma=Theta_Sigma)
  ####################################
  #Sample Precision/Variance Parameter
  ####################################
  # SSE=sum((Obs-process_value-X_M%*%theta[t,])^2)
  # sigma_noise[t]=1/rgamma(1,shape=n_total/2,rate=SSE/2)
  # 
  SSE=sum(unlist(lapply(1:n_unit,FUN=function(i){(Obs[id_index[[i]]]-process_value[id_index[[i]]]-cluster_rd[t-1,cluster_label[id_index[[i]]]]-
                                                    hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
                                                    X_M[id_index[[i]],]%*%theta[t,])^2*kappa[t-1,i]})))
  
  sigma_noise[t]=1/rgamma(1,shape=n_total/2,rate=SSE/2)
  #Sample Kappa
  if(kappa_sample==1){
    for (i in sample(1:n_unit,n_unit))
    {
      kappa_rate=v/2+sum((Obs[id_index[[i]]]-process_value[id_index[[i]]]-cluster_rd[t-1,cluster_label[id_index[[i]]]]-
                            hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
                            X_M[id_index[[i]],]%*%theta[t,])^2)/(2*sigma_noise[t])
      kappa_shape=v/2+length(id_index[[i]])/2
      kappa[t,i]=rgamma(n=1,shape=kappa_shape,rate=kappa_rate)
    }}else{kappa[t,]=kappa[t-1,]}
  
  ##################################
  ##Sample Cluster Random Effect
  ##################################
  for (k in sample(1:length(index_cluster)))
  {
    rd_mean=sum(Obs[index_cluster[[k]]]-process_value[index_cluster[[k]]]-hydro_cluster_rd[t-1,hydro_cluster_label[index_cluster[[k]]]]-
      X_M[index_cluster[[k]],]%*%theta[t,])*cluster_rd_sigma[t-1]/(cluster_rd_sigma[t-1]*length(index_cluster[[k]])+sigma_noise[t])
    rd_var=sigma_noise[t]*cluster_rd_sigma[t-1]/(cluster_rd_sigma[t-1]*length(index_cluster[[k]])+sigma_noise[t])
    cluster_rd[t,k]=rnorm(n=1,mean=rd_mean,sd=sqrt(rd_var))
  }
  ##Sample the Variance for Random Effect
  ##Cluster sigma
  cluster_shape=length(index_cluster)/2+1/2
  cluster_rate=0.5*(sum(cluster_rd[t,]^2)+1)
  cluster_rd_sigma[t]=1/rgamma(n=1,shape=cluster_shape,rate=cluster_rate)
  
  #################################
  ###Sample HydroYear Random Effect
  #################################
  for (k in sample(1:length(hydro_cluster)))
  {
    rd_mean=sum(Obs[hydro_cluster[[k]]]-process_value[hydro_cluster[[k]]]-cluster_rd[t,cluster_label[hydro_cluster[[k]]]]-
                  X_M[hydro_cluster[[k]],]%*%theta[t,])*hydro_cluster_rd_sigma[t-1]/(hydro_cluster_rd_sigma[t-1]*length(hydro_cluster[[k]])+sigma_noise[t])
    rd_var=sigma_noise[t]*hydro_cluster_rd_sigma[t-1]/(hydro_cluster_rd_sigma[t-1]*length(hydro_cluster[[k]])+sigma_noise[t])
    hydro_cluster_rd[t,k]=rnorm(n=1,mean=rd_mean,sd=sqrt(rd_var))
  }
  ##Sample the Variance for Random Effect for HydroYear
  ##Cluster sigma
  cluster_shape=length(hydro_cluster)/2+1/2
  cluster_rate=0.5*(sum(hydro_cluster_rd[t,]^2)+1)
  hydro_cluster_rd_sigma[t]=1/rgamma(n=1,shape=cluster_shape,rate=cluster_rate)
  
  
  ####################################
  #Causal Tau MGP
  ####################################
  if(tau_mgp==1){
    prod=1
    prod_vector=1
    for (k in 2:K)
    {
      prod=prod*tau_delta[t-1,k]
      prod_vector=c(prod_vector,prod)
    }
    
    tau_delta[t,1]=rgamma(1,shape=amu[t-1,1]+K,
                          rate=1+0.5*sum(prod_vector*(tau_0[t,]^2+tau_1[t,]^2)))
    #K>=2
    shuffle_k=sample(2:K,K-1)
    scan=0
    for (k in shuffle_k)
    {
      prod=1
      for (h in 1:(k-1))
      {
        #Update the sampled one
        if (h%in%c(1,shuffle_k[0:scan]))
        {prod=prod*tau_delta[t,h]
        }else{prod=prod*tau_delta[t-1,h]}
        
      }
      prod_vector=prod
      if(k<K){
        for (h in (k+1):K)
        {
          if (h%in%shuffle_k[0:scan])
          {
            prod=prod*tau_delta[t,h]
          }
          else{
            prod=prod*tau_delta[t-1,h]
          }
          prod_vector=c(prod_vector,prod)
        }
      }
      #Truncated Gamma
      if(truncate==1){
        low_u=pgamma(1,shape=amu[t-1,2]+K+1-k,
                     rate=1+0.5*sum(prod_vector*(tau_0[t,k:K]^2+tau_1[t,k:K]^2)))
        u=runif(1,low_u,1)
        tau_delta[t,k]=qgamma(u,shape=amu[t-1,2]+K+1-k,
                              rate=1+0.5*sum(prod_vector*(tau_0[t,k:K]^2+tau_1[t,k:K]^2)))
      }else{tau_delta[t,k]=rgamma(n=1,shape=amu[t-1,2]+K+1-k,
                                  rate=1+0.5*sum(prod_vector*(tau_0[t,k:K]^2+tau_1[t,k:K]^2)))}
      scan=scan+1
    }
    #Reconstruct sigma_tau
    variance=1
    for (k in 1:K)
    {
      
      variance=variance/tau_delta[t,k]
      sigma_tau[t,k]=variance
    }
    
    ###Update amu,ard with MH
    proposal=runif(1,amu[t-1,1]-step_size,amu[t-1,1]+step_size)
    #>0 restriction
    proposal=abs(proposal)
    acceptance_rate=log(tau_delta[t,1])*(proposal-amu[t-1,1])+
      log(gamma(amu[t-1,1]))-log(gamma(proposal))+
      amu[t-1,1]-proposal+log(proposal)-log(amu[t-1,1])
    acceptance_rate=min(exp(acceptance_rate),1)
    if (runif(1)<acceptance_rate)
    {amu[t,1]=proposal}else{amu[t,1]=amu[t-1,1]}
    
    
    proposal=runif(1,amu[t-1,2]-step_size,amu[t-1,2]+step_size)
    #>=2 restriction
    proposal=max(abs(proposal),2)
    acceptance_rate=sum(log(tau_delta[t,2:K]))*(proposal-amu[t-1,2])+
      (log(gamma(amu[t-1,2]))-log(gamma(proposal)))*(K-1)+
      amu[t-1,2]-proposal+log(proposal)-log(amu[t-1,2])
    acceptance_rate=min(exp(acceptance_rate),1)
    if (runif(1)<acceptance_rate)
    {amu[t,2]=proposal}else{amu[t,2]=amu[t-1,2]}
  }else{
    sigma_tau[t]=1/rgamma(n=1,shape=K+1/2,rate=(sum(tau_0[t,]^2+tau_1[t,]^2)+1)/2)
    tau_delta[t,]=tau_delta[t-1,]
    amu[t,]=amu[t-1,]
  }
  
  ####################################
  #Random Effect MGP
  ####################################
  # SSE_Beta=apply(Beta[[t]]-(1-Treatment_Ind)%*%t(tau_0[t,])-Treatment_Ind%*%t(tau_1[t,]),2,
  #                FUN=function(x){sum(x^2)})
  Beta_Error=Beta[[t]]-(1-Treatment_Ind)%*%t(tau_0[t,])-Treatment_Ind%*%t(tau_1[t,])
  #Sample Phi
  if(phi_sample==1){
    phi[[t]]=phi[[1]]
    for (k in sample(1:K,K))
    {
      for (i in sample(1:n_unit,n_unit))
      {
        phi[[t]][i,k]=rgamma(n=1,shape=(v+1)/2,
                             rate=v/2+(Beta_Error[i,k]^2/sigma_K[t-1,k])/2)
        
      }
    }}else{phi[[t]]=phi[[1]]}
  
  
  prod=1
  prod_vector=1
  for (k in 2:K)
  {
    prod=prod*rd_delta[t-1,k]
    prod_vector=c(prod_vector,prod)
  }
  
  rd_delta[t,1]=rgamma(1,shape=ard[t-1,1]+K*n_unit/2,
                       rate=1+0.5*sum(prod_vector*apply(phi[[t]] * Beta_Error^ 2, 2, sum)))
  
  #K>=2
  shuffle_k=sample(2:K,K-1)
  scan=0
  for (k in shuffle_k)
  {
    prod=1
    for (h in 1:(k-1))
    {
      #Update the sampled one
      if (h%in%c(1,shuffle_k[0:scan]))
      {prod=prod*rd_delta[t,h]
      }else{prod=prod*rd_delta[t-1,h]}
      
    }
    prod_vector=prod
    if(k<K){
      for (h in (k+1):K)
      {
        if (h%in%shuffle_k[0:scan])
        {prod=prod*rd_delta[t,h]
        }else{prod=prod*rd_delta[t-1,h]}
        prod_vector=c(prod_vector,prod)
      }
      temp_rate=1+0.5*sum(prod_vector*apply(phi[[t]][,k:K] * Beta_Error[,k:K]^ 2, 2, sum))
      temp_shape=ard[t-1,2]+(K-k+1)*n_unit/2
      if(truncate==1){
        low_u=pgamma(1,shape=temp_shape,rate=temp_rate)
        u=runif(1,low_u,1)
        rd_delta[t,k]=qgamma(u,shape=temp_shape,rate=temp_rate)
      }else{
        rd_delta[t,k]=rgamma(n=1,shape=temp_shape,rate=temp_rate)
      }
      #print(1+0.5*sum(prod_vector*apply(phi[[t]][,k:K] * Beta_Error[,k:K]^ 2, 2, sum)))
    }
    else{
      temp_rate=1+0.5*sum(prod_vector*sum(phi[[t]][,K] * Beta_Error[,K]^ 2))
      temp_shape=ard[t-1,2]+(K-k+1)*n_unit/2
      if(truncate==1){
        low_u=pgamma(1,shape=temp_shape,rate=temp_rate)
        u=runif(1,low_u,1)
        rd_delta[t,k]=qgamma(u,shape=temp_shape,rate=temp_rate)
      }else{
        rd_delta[t,k]=rgamma(n=1,shape=temp_shape,rate=temp_rate)
      }
    }
    scan=scan+1
    
  }
  
  #Reconstruct sigma_K
  variance=1
  for (k in 1:K)
  {
    
    variance=variance/rd_delta[t,k]
    sigma_K[t,k]=variance
  }
  
  ##Update ard with MH
  proposal=runif(1,ard[t-1,1]-step_size,ard[t-1,1]+step_size)
  proposal=abs(proposal)
  acceptance_rate=log(rd_delta[t,1])*(proposal-ard[t-1,1])+
    log(gamma(ard[t-1,1]))-log(gamma(proposal))+
    ard[t-1,1]-proposal+(log(proposal)-log(ard[t-1,1]))*(a_rd_1-1)
  acceptance_rate=min(exp(acceptance_rate),1)
  if (runif(1)<acceptance_rate)
  {ard[t,1]=proposal}else{ard[t,1]=ard[t-1,1]}
  
  
  proposal=runif(1,ard[t-1,2]-step_size,ard[t-1,2]+step_size)
  #  proposal=max(abs(proposal),0)
  proposal=abs(proposal)
  acceptance_rate=sum(log(rd_delta[t,2:K]))*(proposal-ard[t-1,2])+
    (log(gamma(ard[t-1,2]))-log(gamma(proposal)))*(K-1)+
    ard[t-1,2]-proposal+(log(proposal)-log(ard[t-1,2]))*(a_rd_2-1)
  acceptance_rate=min(exp(acceptance_rate),1)
  if (runif(1)<acceptance_rate)
  {ard[t,2]=proposal}else{ard[t,2]=ard[t-1,2]}
  
  ####################################
  if(t%%(simu_time/100)==0){
    print(paste("==",t/simu_time*100,"%=="))}
}

burn_in=seq(simu_time*9/10,simu_time,by=10)
#burn_in=seq(1,simu_time/10,by=10)
pos_sample_size=length(burn_in)
B_Grid_NEW=construct_B(Time_Grid,knots)
#B_Grid_NEW=construct_B_poly(Time_Grid,1:l)
eigen=lapply(1:pos_sample_size,FUN=function(x){B_Grid_NEW%*%(Psi[[burn_in[x]]])})


M_MCMC_Result=list(sigma_K=sigma_K[burn_in,],
                   sigma_noise=sigma_noise[burn_in],
                   tau_0=tau_0[burn_in,],
                   tau_1=tau_1[burn_in,],
                   eigen=eigen,
                   phi=lapply(burn_in,FUN=function(x){phi[[x]]}),
                   sigma_tau=sigma_tau,
                   rd_delta=rd_delta[burn_in,],
                   Beta=lapply(burn_in,FUN=function(x){Beta[[x]]}),
                   Psi=lapply(burn_in,FUN=function(x){Psi[[x]]}),
                   lambda=lambda[burn_in,],
                   theta=theta[burn_in,],
                   cluster_rd=cluster_rd[burn_in,],
                   hydro_cluster_rd=hydro_cluster_rd[burn_in,],
                   cluster_sigma=cluster_rd_sigma[burn_in],
                   knots=knots)

Beta_M=M_MCMC_Result$Beta
Psi_M=M_MCMC_Result$Psi
theta_M=M_MCMC_Result$theta
tau_0_M=M_MCMC_Result$tau_0
tau_1_M=M_MCMC_Result$tau_1
sigma_K=M_MCMC_Result$sigma_K
eigen=M_MCMC_Result$eigen
cluster_rd_M=M_MCMC_Result$cluster_rd
hydro_cluster_rd_M=M_MCMC_Result$hydro_cluster_rd

eigen_mean=apply(simplify2array(eigen),1:2,mean)
eigen_ci_up=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.975)})
eigen_ci_down=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.025)})
# 
order_factor=order(apply(sigma_K,2,mean),decreasing =T)
sort(apply(sigma_K,2,mean),decreasing = T)/sum(apply(sigma_K,2,mean))

# pdf("PC_Compare.pdf",height=6,width=10)
# par(mfrow=c(1,2))
# plot(age_grids,eigen_mean[,order_factor[1]],ylim=range(eigen_mean),type='l',col="red",ylab="Eigen_function",
#      main="Eigen Functions for DSI to Females")
# lines(age_grids,eigen_mean[,order_factor[2]],col="blue")
# #lines(age_grids,eigen_mean[,order_factor[3]],col="green")
# lines(age_grids,eigen_ci_up[,order_factor[1]]+0.1,type='l',lty=2)
# lines(age_grids,eigen_ci_down[,order_factor[1]]-0.1,type='l',lty=2)
# lines(age_grids,eigen_ci_up[,order_factor[2]]+0.2,type='l',lty=2)
# lines(age_grids,eigen_ci_down[,order_factor[2]]-0.2,type='l',lty=2)
# legend("top",legend=c("1st PC 57.66%","2nd PC 25.45%"),lty=1,col=c("red","blue"))
# #lines(age_grids,eigen_ci_up[,order_factor[3]],type='l',lty=2)
# #lines(age_grids,eigen_ci_down[,order_factor[3]],type='l',lty=2)
# 
# #Uncertainty
# # apply(sigma_K[burn_in,],2,mean)
# # apply(sigma_tau[burn_in,],2,mean)
truncate_K=K_m
mediator_effect=lapply(1:pos_sample_size,FUN=function(x){eigen[[x]][,order_factor[c(1,truncate_K)]]%*%(tau_1_M[x,order_factor[c(1,truncate_K)]]-tau_0_M[x,order_factor[c(1,truncate_K)]])})
mediator_effect_mean=apply(simplify2array(mediator_effect),1:2,median)
mediator_effect_down=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.025)})
mediator_effect_up=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.975)})
# 
# # pdf("PC_Compare.pdf",height=6,width=10)
# # par(mfrow=c(1,2))
# # plot(work_grid,eigen_mean[,2],ylim=range(eigen_mean),
# #      type='l',col="red",lwd=1.2,lty=1,ylab="Eigenfunction",
# #      xlab="Time",main="Eigenfunction of Mediators")
# # lines(work_grid,eigen_mean[,3],col="blue",lwd=1.2,lty=1)
# # lines(work_grid,eigen_ci_up[,3],type='l',lty=2)
# # lines(work_grid,eigen_ci_down[,3],type='l',lty=2)
# # lines(work_grid,eigen_ci_up[,2],type='l',lty=2)
# # lines(work_grid,eigen_ci_down[,2],type='l',lty=2)
# # legend("bottomright",legend=c("1st PC","2nd PC"),lty=1,col=c("red","blue"))
# # # 
# # # 
# PC=lapply(1:pos_sample_size,FUN=function(x){Beta_M[[x]][,c(1,2)]})
# unit_PC_Mean=apply(simplify2array(PC),1:2,mean)
# plot(unit_PC_Mean[,1],unit_PC_Mean[,2],xlab="Score of 1st PC",
#      ylab="Score of 2nd PC",col=Treatment_Ind+3,pch=3,cex=0.7,
#      main="Score on 1st and 2nd PC")
# points(mean(unit_PC_Mean[which(Treatment_Ind==0),1]),mean(unit_PC_Mean[which(Treatment_Ind==0),2]))
# points(mean(unit_PC_Mean[which(Treatment_Ind==1),1]),mean(unit_PC_Mean[which(Treatment_Ind==1),2]))
# legend("bottomright",legend=c("Treated","Control"),col=c("blue","green"),
#        pch=3)
#  dev.off()
# points(mean(unit_PC_Mean[which(Treatment_Ind==1),1]),
#        mean(unit_PC_Mean[which(Treatment_Ind==1),2]),cex=3,pch=4)

# pdf(file=paste("M_Effect",K,l,bandwidth,"plot.pdf",sep="_"),width=6,height=4)
# plot(age_grids,mediator_effect_mean-0.1,ylim=range(mediator_effect),type='l',col="red",xlab="Age",ylab="Effect on Mediator",
#      main="Effect on mediator along the life span, DSI_F, females")
# lines(age_grids,mediator_effect_up-0.1,lty=2)
# lines(age_grids,mediator_effect_down-0.1,lty=2)
# abline(h=0)
# dev.off()
