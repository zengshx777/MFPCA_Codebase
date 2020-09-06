
l=13;K=6;Grid.num=40;bandwidth=2;simple=1;simu_time=100000

construct_B<-function(x,knots)
{
  return(cbind(1, x, matrix(unlist(lapply(
    x,
    FUN = function(x) {
      exp(-bandwidth*(x - knots) ^ 2)
    }
  )), nrow = length(x), byrow = T)))
}




time_grids_obs=unique(M_data$Time)
knots_m=quantile(time_grids_obs,seq(0,1,length=l+2))[-c(1,l+2)]

B_list=aggregate(M_data$Time,by=list(M_data$ID),FUN=
                   function(x){construct_B(x,knots_m)})
MB_matrix=B_list$x

Time_Grid=seq(0,1,length=Grid.num)


#Construct Time Grid B
MB_Grid=construct_B(Time_Grid,knots_m)
MB_Grid_Product=t(MB_Grid)%*%MB_Grid


n_unit=length(unique(M_data$ID))
n_total=nrow(M_data)

#Covariate Matrix
X_M=as.matrix(M_data[,c("X1","X2")])
X_M=scale(X_M)
X_M_mean=attr(X_M,"scaled:center")
X_M_sd=attr(X_M,"scaled:scale")
if(simple==1){
  X_M=cbind(rep(1,nrow(X_M)),X_M)
}else{
  X_M=cbind(construct_B(M_data$Time,knots_m),X_M)}
X_Product=t(X_M)%*%X_M

id_index=aggregate(1:nrow(X_M),by=list(M_data$ID),FUN=function(x){x})
id_index=id_index$x


#Outcome
Obs=M_data$Mediator
Obs_list=aggregate(M_data$Mediator,by=
                     list(M_data$ID),
                   function(x){x})
Obs_value=Obs_list$x


Treatment_Ind=aggregate(M_data$Treatment,by=list(M_data$ID),
                        FUN=function(x){mean(x)})$x
N_0=sum(Treatment_Ind==0)
N_1=n_unit-N_0

#Penalty Function
Omega=diag(c(0,0,rep(1,l)))


#MCMC Time
#Parameter for Covariates
theta=matrix(0,ncol=ncol(X_M),nrow=simu_time)
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
sigma_K=matrix(0,ncol=K,nrow=simu_time);sigma_K[1,]=K:1
sigma_theta=10^2
#TE Parameter
tau_0=tau_1=matrix(0,ncol=K,nrow=simu_time)


#Parameter for MGP
#For causal parameters
#tau_delta=matrix(0,ncol=K,nrow=simu_time);tau_delta[1,]=rep(1,K)
#For RD
rd_delta=matrix(0,ncol=K,nrow=simu_time);rd_delta[1,]=seq(1,2,length=K)
sigma_K[1,1]=1
for (k in 2:K)
{
  sigma_K[1,k]=sigma_K[1,k-1]/rd_delta[1,k]
}
sigma_tau=100;
ard=matrix(0,ncol=2,nrow=simu_time);ard[1,]=c(1,3)

phi_ini=matrix(1,ncol=K,nrow=n_unit)
phi=rep(list(NA),simu_time)
phi[[1]]=phi_ini

#Step Size for MH
step_size=1
v=1


##Gibbs Sampler
for (t in 2:simu_time){
  #Fix k
  #Sample Eigen function
  Psi[[t]]=Psi[[t-1]]
  Beta[[t]]=Beta[[t-1]]
  #Shuffle
  for (k in sample(1:K,K)){
    Q_k=Reduce('+',lapply(
      1:n_unit,
      FUN = function(x) {
        t(MB_matrix[[x]]) %*% MB_matrix[[x]] * Beta[[t]][x, k]^2
      }
    ))/sigma_noise[t-1]+lambda[t-1,k]*Omega
    
    
    l_k=Reduce('+',lapply(1:n_unit,FUN=function(i){t(MB_matrix[[i]])%*%(Obs_value[[i]]-X_M[id_index[[i]],]%*%theta[t-1,]-
                                                                          MB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*Beta[[t]][i,k]}))/sigma_noise[t-1]
    #Add a constraint
    C_k=t(cbind(c(1,rep(0,l+1)),Psi[[t]][,-k]))%*%MB_Grid_Product
    
    Q_L=t(chol(Q_k))
    l_tilde=solve(Q_L)%*%l_k
    Psi_0=solve(t(Q_L))%*%(l_tilde+rnorm(l+2))
    
    C_k=t(Psi[[t]][,-k])%*%MB_Grid_Product
    
    C_tilde=solve(Q_L)%*%t(C_k) 
    C_tilde=solve(t(Q_L))%*%C_tilde
    
    Psi[[t]][, k] = Psi_0 - C_tilde%*% 
      solve(C_k %*% C_tilde)%*%C_k%*%Psi_0
    #Norm of f
    temp_norm=as.numeric(t(Psi[[t]][, k])%*%MB_Grid_Product%*%Psi[[t]][, k])
    Psi[[t]][,k]=Psi[[t]][,k]/sqrt(abs(temp_norm))
    Beta[[t]][,k]=Beta[[t]][,k]*sqrt(abs(temp_norm))
    
    ##Sample Lambda
    lambda[t,k]=rgamma(1,shape=(l+1)/2,rate=t(Psi[[t]][,k])%*%Omega%*%Psi[[t]][,k]/2)
    lambda[t,k]=max(10^(-8),lambda[t,k])
  }
  
  ##Sample Principal Score: Beta
  for (i in sample(1:n_unit,n_unit)){
    temp_f=MB_matrix[[i]]%*%Psi[[t]]
    f_sq_sum=apply(temp_f^2,2,sum)
    temp_sigma=1/(f_sq_sum/sigma_noise[t-1]+phi[[t-1]][i,]/sigma_K[t-1,])
    #print(temp_sigma)
    for (k in sample(1:K,K))
    {
      temp_mean=
        sum((Obs_value[[i]]-X_M[id_index[[i]],]%*%theta[t-1,]-
               MB_matrix[[i]]%*%Psi[[t]][,-k]%*%Beta[[t]][i,-k])*temp_f[,k])/sigma_noise[t-1]+
        (tau_0[t-1,k]*(1-Treatment_Ind[i])+tau_1[t-1,k]*Treatment_Ind[i])*phi[[t-1]][i,k]/sigma_K[t-1,k]
      Beta[[t]][i,k]=rnorm(1,mean=temp_mean*temp_sigma[k],sd=sqrt(temp_sigma[k]))
    }
  }
  
  
  
  
  #Sample Tau_0,Tau_1
  tau_0_sigma=1/(apply(phi[[t-1]][Treatment_Ind==0,],2,sum)/sigma_K[t-1,]+1/sigma_tau)
  tau_0_mean=apply(Beta[[t]][Treatment_Ind==0,]*phi[[t-1]][Treatment_Ind==0,],2,sum)/sigma_K[t-1,]
  tau_0[t,]=rnorm(n=K,mean=tau_0_mean*tau_0_sigma,sd=sqrt(tau_0_sigma))
  
  tau_1_sigma=1/(apply(phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]+1/sigma_tau)
  tau_1_mean=apply(Beta[[t]][Treatment_Ind==1,]*phi[[t-1]][Treatment_Ind==1,],2,sum)/sigma_K[t-1,]
  tau_1[t,]=rnorm(n=K,mean=tau_1_mean*tau_1_sigma,sd=sqrt(tau_1_sigma))
  
  #Sample theta
  Theta_Sigma=solve(X_Product/sigma_noise[t-1]+diag(rep(1/sigma_theta,ncol(X_M))))
  
  process_value=unlist(lapply(1:n_unit,FUN=function(x){MB_matrix[[x]]%*%Psi[[t]]%*%Beta[[t]][x,]}))
  Theta_Mean=Theta_Sigma%*%(t(X_M)%*%(Obs-process_value)/sigma_noise[t-1])
  
  theta[t,]=mvrnorm(n=1,mu=Theta_Mean,Sigma=Theta_Sigma)
  
  #Sample Precision/Variance Parameter
  SSE=sum((Obs-process_value-X_M%*%theta[t,])^2)
  sigma_noise[t]=1/rgamma(1,shape=n_total/2,rate=SSE/2)

  Beta_Error=Beta[[t]]-(1-Treatment_Ind)%*%t(tau_0[t,])-Treatment_Ind%*%t(tau_1[t,])
  #Sample Phi
  phi[[t]]=phi[[1]]
  for (k in sample(1:K,K))
  {
    for (i in sample(1:n_unit,n_unit))
    {
      phi[[t]][i,k]=rgamma(n=1,shape=(v+1)/2,
                           rate=v/2+(Beta_Error[i,k]^2/sigma_K[t-1,k])/2)
      
    }
  }
  
  
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
      
      rd_delta[t,k]=rgamma(1,shape=ard[t-1,2]+(K-k+1)*n_unit/2,
                           rate=1+0.5*sum(prod_vector*apply(phi[[t]][,k:K] * Beta_Error[,k:K]^ 2, 2, sum)))
      #print(1+0.5*sum(prod_vector*apply(phi[[t]][,k:K] * Beta_Error[,k:K]^ 2, 2, sum)))
    }
    else{
      rd_delta[t,k]=rgamma(1,shape=ard[t-1,2]+(K-k+1)*n_unit/2,
                           rate=1+0.5*sum(prod_vector*sum(phi[[t]][,K] * Beta_Error[,K]^ 2)))
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
  ##Update amu,ard with MH
  proposal=runif(1,ard[t-1,1]-step_size,ard[t-1,1]+step_size)
  proposal=abs(proposal)
  acceptance_rate=log(rd_delta[t,1])*(proposal-ard[t-1,1])+
    log(gamma(ard[t-1,1]))-log(gamma(proposal))+
    ard[t-1,1]-proposal+log(proposal)-log(ard[t-1,1])
  acceptance_rate=min(exp(acceptance_rate),1)
  if (runif(1)<acceptance_rate)
  {ard[t,1]=proposal}else{ard[t,1]=ard[t-1,1]}
  
  
  proposal=runif(1,ard[t-1,2]-step_size,ard[t-1,2]+step_size)
  proposal=max(abs(proposal),2)
  acceptance_rate=sum(log(rd_delta[t,2:K]))*(proposal-ard[t-1,2])+
    (log(gamma(ard[t-1,2]))-log(gamma(proposal)))*(K-1)+
    ard[t-1,2]-proposal+log(proposal)-log(ard[t-1,2])
  acceptance_rate=min(exp(acceptance_rate),1)
  if (runif(1)<acceptance_rate)
  {ard[t,2]=proposal}else{ard[t,2]=ard[t-1,2]}
  
  
  if(t%%1000==0){
    print(paste("==",t,"=="))}
  
  
}

# M_MCMC_Result=list(sigma_K=sigma_K,sigma_noise=sigma_noise,
#                    tau_0=tau_0,tau_1=tau_1,Beta=Beta,
#                    Psi=Psi,lambda=lambda,theta=theta,
#                    knots_m=knots_m)

Beta_M=Beta
Psi_M=Psi
theta_M=theta
tau_0_M=tau_0
tau_1_M=tau_1

pos_sample_size=simu_time/10
burn_in=(simu_time*9/10+1):simu_time
B_Grid_NEW=construct_B(seq(0,1,length=200),knots_m)
eigen=lapply(1:pos_sample_size,FUN=function(x){B_Grid_NEW%*%(Psi_M[[burn_in[x]]])})
eigen_mean=apply(simplify2array(eigen),1:2,mean)
eigen_ci_up=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.975)})
eigen_ci_down=apply(simplify2array(eigen),1:2,FUN=function(x){quantile(x,0.025)})


# plot(eigen_mean[,3],ylim=range(eigen_mean),type='l',col="red")
# lines(eigen_mean[,4],col="blue")
# lines(eigen_ci_up[,3],type='l',lty=2)
# lines(eigen_ci_down[,3],type='l',lty=2)
# lines(eigen_ci_up[,4],type='l',lty=2)
# lines(eigen_ci_down[,4],type='l',lty=2)
# lines(eigen_mean[,1],col="blue")
# lines(eigen_ci_up[,1],type='l',lty=2)
# lines(eigen_ci_down[,1],type='l',lty=2)
# lines(eigen_mean[,2],col="blue")
# lines(eigen_ci_up[,2],type='l',lty=2)
# lines(eigen_ci_down[,2],type='l',lty=2)

#Uncertainty
# #apply(sigma_tau[burn_in,],2,mean)
# apply(sigma_K[burn_in,],2,mean)
# apply(rd_delta[burn_in,],2,mean)

mediator_effect=lapply(1:pos_sample_size,FUN=function(x){eigen[[x]][,]%*%(tau_1_M[burn_in[x],]-tau_0_M[burn_in[x],])})
mediator_effect_mean=apply(simplify2array(mediator_effect),1:2,mean)
mediator_effect_down=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.015)})
mediator_effect_up=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.965)})

mediator_effect=lapply(1:pos_sample_size,FUN=function(x){eigen[[x]][,]%*%(tau_1_M[burn_in[x],]-tau_0_M[burn_in[x],])})
mediator_effect_mean=apply(simplify2array(mediator_effect),1:2,mean)
mediator_effect_down=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.025)})
mediator_effect_up=apply(simplify2array(mediator_effect),1:2,function(x){quantile(x,0.975)})

pdf("Effect_on_M_No.pdf",width=8,height=6)
plot(work_grid,mediator_effect_mean,ylim=range(mediator_effect),xlab="Time",ylab="Effect on Mediator",
     type='l',lwd=2,lty=5)
lines(work_grid,mediator_effect_up,lty=2)
lines(work_grid,mediator_effect_down,lty=2)
lines(work_grid,m_1-m_0,lwd=2)
abline(h=0)
legend("top",legend=c("True Value","Posterior Mean","Credible Interval"),lty=c(1,5,2),lwd=c(2,2,1))
dev.off()
# ##Recover
# #Impute M_Process Value
t_rd=sample((simu_time*9/10+1):simu_time,1)
#Use the Parameter from M estimation
m_process = lapply((simu_time*9/10+1):simu_time,
      FUN=function(y){
      MB_Grid %*% Psi_M[[y]] %*% t(Beta_M[[y]])
    })
m_process_mean=apply(simplify2array(m_process),1:2,mean)

# 
mean_process=apply(cbind(rep(1,nrow(X_M)),X_M[,c("X1","X2")])%*% t(theta_M),1,mean)
 random_id=sample(1:200,1)
# random_id=151
# pdf("M_Process_Recovery.pdf",width=8,height=6)
plot(M_data[id_index[[random_id]],"Time"],M_data[id_index[[random_id]],"True_Mediator"]-mean_process[id_index[[random_id]]],type='o',
#       -(-M_data$X1[id_index[[random_id]]]+0.5*M_data$X2[id_index[[random_id]]]),col="blue",type='o',
     xlab="Time",ylab="Mediator Value",xlim=c(0,1),ylim=range(M_data$Mediator))
lines(Time_Grid,m_process_mean[,random_id],lty=2)
#lines(M_data[id_index[[random_id]],"Time"],m_process[id_index[[random_id]]],col="red")
legend("top",legend=c("True","Posterior mean"),lty=c(1,2))
#dev.off()
