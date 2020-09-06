eta_M=M_MCMC_Result$Beta
Psi_M=M_MCMC_Result$Psi
theta_M=M_MCMC_Result$theta
tau_0_M=M_MCMC_Result$tau_0
tau_1_M=M_MCMC_Result$tau_1
sigma_K=M_MCMC_Result$sigma_K
cluster_rd_M=M_MCMC_Result$cluster_rd
hydro_cluster_rd_M=M_MCMC_Result$hydro_cluster_rd

#library(expm)
# construct_B<-function(x,knots)
# {
#   #Correlation
#   Omega_Z0=matrix(unlist(lapply(
#     knots,FUN=function(x){abs((x-knots)^3)}
#   )),nrow=length(knots),byrow=T)
#   diag(Omega_Z0)=1
#   return(cbind(1, x, matrix(unlist(lapply(
#     x,
#     FUN = function(x) {
#       abs((x - knots) ^ 3)
#     }
#   )), nrow = length(x), byrow = T)%*%sqrtm(solve(Omega_Z0)) ))
# }


#New Knots, maybe different with the mediator
knots=quantile(Complete.data$t,seq(0,1,length=l+2))[-c(1,l+2)]

B_list=aggregate(Complete.data$t,by=list(Complete.data$sname),FUN=
                   function(x){construct_B(x,knots)})
YB_matrix=B_list$x
# Time_Grid=seq(0,1,length=Grid.num)

#Construct Time Grid B
YB_Grid=construct_B(Time_Grid,knots)
YB_Grid_Product=t(YB_Grid)%*%YB_Grid


#Covariate Matrix
# X_Y = model.matrix(as.formula(paste(
#   "~", paste(covariate.index, collapse = "+")
# )), Complete.data)

X_Y = model.matrix(as.formula(paste(
  "~", paste(covariate.index, collapse = "+")
)), Complete.data)
X_Y=X_Y[,-1]
#Centralize
X_Y=scale(X_Y)
X_Y_mean=attr(X_Y,"scaled:center")
X_Y_sd=attr(X_Y,"scaled:scale")
if(simple==1){X_Y=cbind(rep(1,nrow(X_Y)),X_Y)
}else{
  X_Y=cbind(construct_B(Complete.data$t,knots),X_Y)}
#AUG_X=cbind(X_Y,Complete.data[,mediation.index[m]])

id_index=aggregate(1:nrow(X_Y),by=list(Complete.data$sname),FUN=function(x){x})
id_index=id_index$x

# ##Cluster ID
# # cluster_label_value=unique(Complete.data$grp)
# # cluster_label=1:length(cluster_label_value)
# index_cluster=aggregate(1:nrow(Complete.data),by=list(Complete.data$grp),FUN=function(x){x})
# #12 cluster in total
# ##index label for each cluster
# cluster_label=unlist(lapply(Complete.data$grp,FUN=function(x){which(index_cluster$Group.1==x)}))
# ##The index set for each cluster
# index_cluster=index_cluster$x



#Outcome
Obs=Complete.data$gc
Obs_list=aggregate(Complete.data$gc,by=
                     list(Complete.data$sname),
                   function(x){x})
Obs_value=Obs_list$x

#Treatment Indicator
Treatment_Ind=aggregate(Complete.data[,adverse.index[adv_id]],by=list(Complete.data$sname),
                        FUN=function(x){if(any(x)){return(1)}else{return(0)}})$x
N_0=sum(Treatment_Ind==0)
N_1=n_unit-N_0

#Penalty Function
Omega=diag(c(0,0,rep(1,l)))

#MCMC Time
#simu_time=20000
#Parameter for Covariates
#Parameter for Covariates
theta=matrix(0,ncol=ncol(X_Y)+1,nrow=simu_time)
#Principal Score
Beta_ini=matrix(runif(K*n_unit,-0.5,0.5),ncol=K,nrow=n_unit)
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

sigma_theta=10^2

#TE Parameter
tau_0=tau_1=matrix(0,ncol=K,nrow=simu_time)

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


##Cluster Random Effect
cluster_rd=matrix(0,ncol=length(index_cluster),nrow=simu_time)
cluster_rd_sigma=rep(1,simu_time)

hydro_cluster_rd=matrix(0,ncol=length(hydro_cluster),nrow=simu_time)
hydro_cluster_rd_sigma=rep(1,simu_time)

#rd_xi=matrix(0,ncol=n_unit,nrow=simu_time)
#sigma_tau=matrix(0,ncol=K,nrow=simu_time);sigma_tau[1,]=rep(1,K)
#sigma_tau=100
#Hyper
amu=ard=matrix(0,ncol=2,nrow=simu_time);amu[1,]=c(1,2)
ard=matrix(0,ncol=2,nrow=simu_time);ard[1,]=c(1,2)

phi_ini=matrix(1,ncol=K,nrow=n_unit)
phi=rep(list(NA),simu_time)
phi[[1]]=phi_ini

#Step Size for MH
step_size=2
v=5
#HP for MGP
a_rd_1=2
a_rd_2=3

#Impute M_Process Value and Multiply with the spline basis functions
m_estimate=lapply(1:pos_sample_size,FUN=function(t_rd){
  #Use the Parameter from M estimation
  if(simple!=1){
    m_process = unlist(lapply(
      1:n_unit,
      FUN = function(x) {
        MB_matrix[[x]] %*% Psi_M[[t_rd]] %*% Beta_M[[t_rd]][x, ]+cluster_rd_M[t_rd,cluster_label[id_index[[x]]]]+
          +hydro_cluster_rd_M[t_rd,hydro_cluster_label[id_index[[x]]]]
      }
    )) + cbind(X_Y[,1:(ncol(X_M)-length(covariate.index.dsi))],X_M[,-1:(ncol(X_M)-length(covariate.index.dsi))])%*% theta_M[t_rd, ]
    
  }
  else{
    m_process = unlist(lapply(
      1:n_unit,
      FUN = function(x) {
        MB_matrix[[x]] %*% Psi_M[[t_rd]] %*% Beta_M[[t_rd]][x, ]+cluster_rd_M[t_rd,cluster_label[id_index[[x]]]]+
          hydro_cluster_rd_M[t_rd,hydro_cluster_label[id_index[[x]]]]
      }
    )) + cbind(rep(1,nrow(X_Y)),X_M[,-1])%*% theta_M[t_rd, ]
    
  }
  return(m_process)})
m_plug_in=Reduce("+",m_estimate)/length(burn_in)


#AUG_X=cbind(X_Y,m_plug_in)
#Matrix for M-process
#AUG_X=cbind(X_Y,Complete.data[,mediation.index[m]])
AUG_X=cbind(X_Y,m_plug_in)
X_Product=t(AUG_X)%*%AUG_X

##Gibbs Sampler
for (t in 2:simu_time){
  
  #Impute M_Process Value
  # t_rd=sample((simu_time/2+1):simu_time,1)
  # m_process = unlist(lapply(
  #   1:n_unit,
  #   FUN = function(x) {
  #     B_matrix[[x]] %*% Psi_M[[t_rd]] %*% Beta_M[[t_rd]][x, ]
  #   }
  # )) + X_M %*% theta_M[t_rd, ]
  #AUG_X=cbind(X_Y,DSI_F)
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
        t(YB_matrix[[x]]) %*% YB_matrix[[x]] * Beta[[t]][x, k]^2
      }
    ))/sigma_noise[t-1]+lambda[t-1,k]*Omega
    
    
    l_k=Reduce('+',lapply(1:n_unit,FUN=function(i){t(YB_matrix[[i]])%*%(Obs_value[[i]]-AUG_X[id_index[[i]],]%*%theta[t-1,]-
                                                                          cluster_rd[t-1,cluster_label[id_index[[i]]]]-
                                                                          hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
                                                                          YB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*Beta[[t]][i,k]}))/sigma_noise[t-1]
    #Add a constraint
    C_k=t(cbind(c(1,rep(0,l+1)),Psi[[t]][,-k]))%*%YB_Grid_Product
    
    Q_L=t(chol(Q_k))
    l_tilde=solve(Q_L)%*%l_k
    Psi_0=solve(t(Q_L))%*%(l_tilde+rnorm(l+2))
    
    C_k=t(Psi[[t]][,-k])%*%YB_Grid_Product
    
    C_tilde=solve(Q_L)%*%t(C_k) 
    C_tilde=solve(t(Q_L))%*%C_tilde
    
    Psi[[t]][, k] = Psi_0 - C_tilde%*% 
      solve(C_k %*% C_tilde)%*%C_k%*%Psi_0
    
    temp_norm=abs(as.numeric(t(Psi[[t]][, k])%*%YB_Grid_Product%*%Psi[[t]][, k]))
    Psi[[t]][,k]=Psi[[t]][,k]/sqrt(temp_norm)
    Beta[[t]][,k]=Beta[[t]][,k]*sqrt(temp_norm)
    
    ##Sample Lambda
    lambda[t,k]=rgamma(1,shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    if(shrink==1)
    {
      lambda[t,k]=max(sigma_K[t-1,k],lambda[t,k])
    }else{
      lambda[t,k]=max(10^(-8),lambda[t,k])}
  }
  
  ####################################
  ##Sample Principal Score: Beta
  ####################################
  for (i in sample(1:n_unit,n_unit)){
    temp_f=YB_matrix[[i]]%*%Psi[[t]]
    f_sq_sum=apply(temp_f^2,2,sum)
    temp_sigma=1/(f_sq_sum/sigma_noise[t-1]+phi[[t-1]][i,]/sigma_K[t-1,])
    for (k in sample(1:K,K))
    {
      temp_mean=
        sum((Obs_value[[i]]-AUG_X[id_index[[i]],]%*%theta[t-1,]-cluster_rd[t-1,cluster_label[id_index[[i]]]]-
               hydro_cluster_rd[t-1,hydro_cluster_label[id_index[[i]]]]-
               YB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*temp_f[,k])/sigma_noise[t-1]+
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
  Theta_Sigma = solve(X_Product / sigma_noise[t - 1] + diag(rep(1 / sigma_theta, ncol(AUG_X))))
  
  process_value = unlist(lapply(
    1:n_unit,
    FUN = function(x) {
      YB_matrix[[x]] %*% Psi[[t]] %*% Beta[[t]][x, ]
    }
  ))
  Theta_Mean = Theta_Sigma %*% (t(AUG_X) %*% (Obs - process_value -
                                                unlist(lapply(
                                                  1:n_unit,
                                                  FUN = function(i) {
                                                    cluster_rd[t - 1, cluster_label[id_index[[i]]]] + hydro_cluster_rd[t - 1, hydro_cluster_label[id_index[[i]]]]
                                                  }
                                                ))) / sigma_noise[t - 1])
  
  theta[t, ] = mvrnorm(n = 1, mu = Theta_Mean, Sigma = Theta_Sigma)
  
  ####################################
  #Sample Precision/Variance Parameter
  ####################################
  SSE = sum((Obs - process_value - AUG_X %*% theta[t, ] -
               unlist(lapply(
                 1:n_unit,
                 FUN = function(i) {
                   cluster_rd[t - 1, cluster_label[id_index[[i]]]] + hydro_cluster_rd[t - 1, hydro_cluster_label[id_index[[i]]]]
                 }
               ))) ^ 2)
  sigma_noise[t] = 1 / rgamma(1, shape = n_total / 2, rate = SSE / 2)
  
  
  ##################################
  ##Sample Cluster Random Effect
  ##################################
  for (k in sample(1:length(index_cluster)))
  {
    rd_mean=sum(Obs[index_cluster[[k]]]-process_value[index_cluster[[k]]]-hydro_cluster_rd[t-1,hydro_cluster_label[index_cluster[[k]]]]-
                  AUG_X[index_cluster[[k]],]%*%theta[t,])*cluster_rd_sigma[t-1]/(cluster_rd_sigma[t-1]*length(index_cluster[[k]])+sigma_noise[t])
    rd_var=sigma_noise[t]*cluster_rd_sigma[t-1]/(cluster_rd_sigma[t-1]*length(index_cluster[[k]])+sigma_noise[t])
    cluster_rd[t,k]=rnorm(n=1,mean=rd_mean,sd=sqrt(rd_var))
  }
  ##Sample the Variance for Random Effect
  ##Cluster sigma
  cluster_shape=length(index_cluster)/2+1/2
  cluster_rate=0.5*(sum(cluster_rd[t,]^2)+1)
  cluster_rd_sigma[t]=1/rgamma(n=1,shape=cluster_shape,rate=cluster_rate)
  
  ##################################
  ##Sample Hydro Cluster Random Effect
  ##################################
  for (k in sample(1:length(hydro_cluster)))
  {
    rd_mean=sum(Obs[hydro_cluster[[k]]]-process_value[hydro_cluster[[k]]]-cluster_rd[t,cluster_label[hydro_cluster[[k]]]]-
                  AUG_X[hydro_cluster[[k]],]%*%theta[t,])*cluster_rd_sigma[t-1]/(cluster_rd_sigma[t-1]*length(hydro_cluster[[k]])+sigma_noise[t])
    rd_var=sigma_noise[t]*hydro_cluster_rd_sigma[t-1]/(hydro_cluster_rd_sigma[t-1]*length(hydro_cluster[[k]])+sigma_noise[t])
    hydro_cluster_rd[t,k]=rnorm(n=1,mean=rd_mean,sd=sqrt(rd_var))
  }
  ##Sample the Variance for Random Effect
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
pos_sample_size=length(burn_in)
YB_Grid_NEW=construct_B(Time_Grid,knots)
eigen=lapply(1:pos_sample_size,FUN=function(x){YB_Grid_NEW%*%(Psi[[burn_in[x]]])})

Y_MCMC_Result=list(sigma_K=sigma_K[burn_in,],
                   sigma_noise=sigma_noise[burn_in],
                   phi=lapply(burn_in,FUN=function(x){phi[[x]]}),
                   sigma_tau=sigma_tau,
                   tau_0=tau_0[burn_in,],
                   tau_1=tau_1[burn_in,],
                   eigen=eigen,
                   rd_delta=rd_delta[burn_in,],
                   Beta=lapply(burn_in,FUN=function(x){Beta[[x]]}),
                   Psi=lapply(burn_in,FUN=function(x){Psi[[x]]}),
                   lambda=lambda[burn_in,],
                   theta=theta[burn_in,],
                   cluster_rd=cluster_rd[burn_in,],
                   hydro_cluster_rd=hydro_cluster_rd[burn_in,],
                   hydro_cluster_sigma=hydro_cluster_rd_sigma[burn_in],
                   cluster_sigma=cluster_rd_sigma[burn_in],
                   knots=knots)

#save(Y_MCMC_Result,file=paste(mediation.index[m],adverse.index[adv_id],f,"Y0428.RData",sep="_"))
Beta=Y_MCMC_Result$Beta
Psi=Y_MCMC_Result$Psi
theta=Y_MCMC_Result$theta
tau_0=Y_MCMC_Result$tau_0
tau_1=Y_MCMC_Result$tau_1
sigma_K=Y_MCMC_Result$sigma_K

eigen_mean=apply(simplify2array(eigen),1:2,mean)
eigen_ci_up=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.975)})
eigen_ci_down=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.025)})

order_factor=order(apply(sigma_K,2,mean),decreasing =T)
# pdf("PC_Compare_y.pdf",height=6,width=10)
# par(mfrow=c(1,2))
# plot(age_grids,eigen_mean[,order_factor[1]],ylim=range(eigen_mean),type='l',col="red",ylab="Eigen_function",
#      main="Eigen Functions for Logged GC",xlab="age")
# lines(age_grids,eigen_mean[,order_factor[2]],col="blue")
# #lines(age_grids,eigen_mean[,order_factor[3]],col="green")
# lines(age_grids,eigen_ci_up[,order_factor[1]]+0.02,type='l',lty=2)
# lines(age_grids,eigen_ci_down[,order_factor[1]]-0.02,type='l',lty=2)
# lines(age_grids,eigen_ci_up[,order_factor[2]]+0.03,type='l',lty=2)
# lines(age_grids,eigen_ci_down[,order_factor[2]]-0.03,type='l',lty=2)
# legend("top",legend=c("1st PC 54.37%","2nd PC 24.42%"),lty=1,col=c("red","blue"))
# PC=lapply(1:pos_sample_size,FUN=function(x){Beta[[x]][,c(1,2)]})
# unit_PC_Mean=apply(simplify2array(PC),1:2,mean)
# plot(unit_PC_Mean[,1],unit_PC_Mean[,2],xlab="Score of 1st PC",
#      ylab="Score of 2nd PC",col=Treatment_Ind+3,pch=3,cex=0.7,
#      main="Score on 1st and 2nd PC")
# points(mean(unit_PC_Mean[which(Treatment_Ind==0),1]),mean(unit_PC_Mean[which(Treatment_Ind==0),2]))
# points(mean(unit_PC_Mean[which(Treatment_Ind==1),1]),mean(unit_PC_Mean[which(Treatment_Ind==1),2]))
# legend("bottomright",legend=c("Treated","Control"),col=c("blue","green"),
#        pch=3)
# dev.off()
#Uncertainty
#apply(sigma_K[burn_in,],2,mean)

##Direct Effect
direct_process=lapply(1:pos_sample_size,FUN=function(x){eigen[[x]]%*%(tau_1[x,]-tau_0[x,])})
direct_process_mean=apply(simplify2array(direct_process),1:2,mean)
direct_process_down=apply(simplify2array(direct_process),1:2,function(x){quantile(x,0.025)})
direct_process_up=apply(simplify2array(direct_process),1:2,function(x){quantile(x,0.975)})
# 
# plot(direct_process_mean,ylim=range(direct_process),type='l',col="red")
# lines(direct_process_up,lty=2)
# lines(direct_process_down,lty=2)
# abline(h=0)

##Mediation Effect
gamma=theta[,ncol(theta)]
indirect_process=lapply(1:pos_sample_size,FUN=function(x){
  gamma[x]*mediator_effect[[x]]})
indirect_process_mean=apply(simplify2array(indirect_process),1:2,mean)
indirect_process_down=apply(simplify2array(indirect_process),1:2,function(x){quantile(x,0.025)})
indirect_process_up=apply(simplify2array(indirect_process),1:2,function(x){quantile(x,0.975)})

# indirect_process_mean=mean(gamma)*(mediator_effect_mean-0.1)
# indirect_process_up=mean(gamma)*(mediator_effect_up-0.1)
# indirect_process_down=mean(gamma)*(mediator_effect_down-0.1)

total_effect=lapply(1:pos_sample_size,FUN=function(x){
  indirect_process[[x]]+direct_process[[x]]})
total_effect_mean=apply(simplify2array(total_effect),1:2,mean)
total_effect_down=apply(simplify2array(total_effect),1:2,function(x){quantile(x,0.025)})
total_effect_up=apply(simplify2array(total_effect),1:2,function(x){quantile(x,0.975)})

