library(MASS)
gaussprocess <- function(grids,K = function(s, t) {min(s, t)},
                         start = NULL) {
  t <- grids
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = length(grids)), Sigma = Sigma)
  if (!is.null(start)) {
    path <- path - path[1] + start  # Must always start at "start"
  }
  
  return(data.frame("t" = t, "xt" = path))
}

#Ground Truth
#Potential value of Mediation Process
work_grid=seq(0,1,length=200)
m_0=2*work_grid+sin(work_grid*2*pi)
m_1=4*work_grid+2*sin(2*work_grid*pi)+0.2
# plot(work_grid,m_0,type='l',col="blue",ylim=range(m_0,m_1))
# lines(work_grid,m_1,col="red")
# 
# plot(work_grid,m_1-m_0,type='l',main="Effect on Mediator")
# 
y_1_1=m_1+0.3*work_grid^2+5*work_grid+2*cos(work_grid*2*pi)
y_0_1=m_1+cos(work_grid*2*pi)+0.1*work_grid^2+2*work_grid
y_0_0=m_0+cos(work_grid*2*pi)+0.1*work_grid^2+2*work_grid

# plot(work_grid,y_0_0,type='l',ylim=range(y_0_0,y_0_1,y_1_1),col="blue")
# lines(work_grid,y_0_1,col="green")
# lines(work_grid,y_1_1,col="red")

pdf("True_Value.pdf",width=8,height=6)
plot(work_grid,y_1_1-y_0_0,type='l',ylim=c(0,6),
     xlab="Time",ylab="Total/mediation Effect",lty=1,lwd=2)
lines(work_grid,y_0_1-y_0_0,lty=6,lwd=2)
legend("top",legend=c("Total Effect","Mediation Effect"),lty=c(1,6))
dev.off()
Simulation_Data<-function(nsample_m_mean=20,
                          nsample_y_mean=20,
                          nsample_size=100,
                          same_grids=T,prob_treated=0.5,
                          #Noise for Underlying Process
                          noise_m=0.2,
                          noise_y=0.2,
                          #Noise for Sampling
                          noise_sample=0.5,
                          #Correlation Between T and M
                          p_t_m=0,
                          #Correlation Between T and Y
                          p_t_y=0,
                          #Correlation Between M and Y
                          p_m_y=0
){
  #Sample the observation points
  
  #Same Grid for Mediation and Outcome Process
  nsample_m=rpois(nsample_size,lambda=nsample_m_mean)
  
  if(same_grids==T){
    nsample_y=nsample_m
  }else{
    nsample_y=rpois(nsample_size,lambda=nsample_y_mean)
  }
  
  #Treatment Indicator
  #Equal Probability of Being Treated Or Control
  #Error Term for T,Y,M
  Sigma_Confounding=matrix(c(1,p_t_m,p_t_y,
                             p_t_m,1,p_m_y,
                             p_t_y,p_m_y,1),3,3)

  #Use an Error Term
  Confounding_Error=mvrnorm(n=nsample_size,mu=c(0,0,0),
                            Sigma=Sigma_Confounding)

  treated=Confounding_Error[,1]<qnorm(prob_treated)
  
  
  M_data=NULL
  Y_data=NULL

  for (i in 1:nsample_size)
  {

    
    # #Same Grids
    if(same_grids){
      grids_union=sort(runif(nsample_m[i]))
      while(length(grids_union)!=unique(length(grids_union)))
      {
        grids_union=sort(runif(nsample_m[i]))
      }
      grids_m=grids_y=1:length(grids_union)
      ## Different Grids
    }else{
      # #Potential Grids for Mediator and Outcome Process
      grids_union=sort(runif(nsample_m[i]+nsample_y[i]))
      while(length(grids_union)!=unique(length(grids_union)))
      {
        grids_union=sort(runif(nsample_m[i]+nsample_y[i]))
      }
      #Grids for Mediator Process
      grids_m=sort(sample(length(grids_union),nsample_m[i],replace=F))
      #Grids For Outcome Process
      grids_y=sort(sample(length(grids_union),nsample_y[i],replace=F))
    }
    
    #Underlying Potential Values for Mediator
    mtrue_0=2*grids_union+sin(grids_union*2*pi)
    mtrue_1=4*grids_union+2*sin(grids_union*2*pi)+0.2
    #Covariates
    X=matrix(rnorm(n=length(grids_union)*3,sd=1),ncol=3,nrow=length(grids_union))
    #Obs Mediators
    m_true=mtrue_0*(1-treated[i])+mtrue_1*treated[i]+X[,c(1,2)]%*%c(-1,0.5)+
      gaussprocess(grids_union,
                   K=function(s, t) {noise_m*exp(-8*(s-t)^2)})$xt+Confounding_Error[i,2]*noise_m
    
    
    # plot(grids_union[grids_m],mtrue_0[grids_m],ylim=range(mtrue_0,mtrue_1),col="red",type='l')
    # lines(grids_union[grids_m],mtrue_1[grids_m],col="blue")
    # points(grids_union[grids_m],m_obs[grids_m])
    
    #Underlying Potential Values for Outcome Process
    ytrue_0_0=mtrue_0+cos(grids_union*2*pi)+0.1*grids_union^2+2*grids_union
    ytrue_0_1=mtrue_1+cos(grids_union*2*pi)+0.1*grids_union^2+2*grids_union
    ytrue_1_1=mtrue_1+2*cos(grids_union*2*pi)+
      0.3*grids_union^2+5*grids_union
    
    #First, evaluate the mediator process on y-grids
    #y_true=m_true+cos(grids_union*2*pi)+0.2*grids_union^2+3*grids_union+cos(grids_union*2*pi)
    
    y_true=m_true+cos(grids_union*2*pi)+0.1*grids_union^2+2*grids_union+
      treated[i]*(0.2*grids_union^2+3*grids_union+cos(grids_union*2*pi))+X[,c(2,3)]%*%c(-0.5,1)+
      gaussprocess(grids_union,
                   K=function(s, t) {noise_y*exp(-8*(s-t)^2)})$xt+Confounding_Error[i,3]*noise_y
    
    # plot(grids_union,ytrue_0_0,ylim=range(ytrue_0_0,ytrue_0_1,ytrue_1_1),type='l',col="red")
    # lines(grids_union,ytrue_0_1,col="green")
    # lines(grids_union,ytrue_1_1,col="blue")
    # points(grids_union,yobs)
    
    m_obs=m_true+rnorm(length(grids_union),0,noise_sample)
    y_obs=y_true+rnorm(length(grids_union),0,noise_sample)
    
    M_data=rbind(M_data,cbind(rep(i,nsample_m[i]),rep(treated[i],nsample_m[i]),
                              grids_union[grids_m],m_obs[grids_m],m_true[grids_m],X[grids_m,]))

    
    Y_data=rbind(Y_data,cbind(rep(i,nsample_y[i]),rep(treated[i],nsample_y[i]),
                              grids_union[grids_y],y_obs[grids_y],y_true[grids_y],m_true[grids_y],X[grids_y,]))
  }
  
  colnames(M_data)=c("ID","Treatment","Time","Mediator","True_Mediator","X1","X2","X3")
  colnames(Y_data)=c("ID","Treatment","Time","Outcome","True_Outcome","True_Mediator","X1","X2","X3")
  M_data<-as.data.frame(M_data)
  Y_data<-as.data.frame(Y_data)
  
  treated_id=which(treated==1)
  control_id=which(treated==0)
  
  return(list(M_data=M_data,
              Y_data=Y_data,
              treated_id=treated_id,
              control_id=control_id))
}

