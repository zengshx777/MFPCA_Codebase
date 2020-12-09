rm(list=ls())


args=commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#Illustration
source("Simulation_Function.R")
library(mgcv)
library(gee)

# Simulated dataset
for (i in 1:1000){
Simu_Wrap=Simulation_Data(nsample_m_mean=sparsity.level,
                          nsample_y_mean=sparsity.level,
                          nsample_size=sample.size,
                          same_grids=T,prob_treated=0.5,
                          #Noise for Underlying Process
                          noise_m=1,
                          noise_y=1,
                          #Noise for Sampling
                          noise_sample=2,
                          p_t_m = 0,
                          p_t_y = 0,
                          p_m_y = 0)
M_data=Simu_Wrap$M_data
Y_data=Simu_Wrap$Y_data
source("Simu_Mediator.R")
source("Simu_Sampling_Outcome.R")


true_direct_effect_process = mean(y_1_1-y_0_1)
true_indirect_effect_process = mean(y_0_1-y_0_0)
true_total_effect_process = mean(y_1_1-y_0_0)
gamma_true =1 
direct_effect_process=apply(simplify2array(direct_process),1:2,mean)
indirect_effect_process=apply(simplify2array(indirect_process),1:2,mean)
total_effect_process=apply(simplify2array(total_effect),1:2,mean)
gamma_est = mean(gamma)


## RANDOM EFFECT MODEL
m_model_re=gam(Mediator~X1+X2+X3+Treatment+s(Time)+s(ID,bs="re") ,data=M_data)
Y_data$Mediator = M_data$Mediator
y_model_re=gam(Outcome~X1+X2+X3+Mediator+Treatment+s(Time)+s(ID,bs="re") ,data=Y_data)



## GEE MODEL
gee_model_1=gee(Mediator~X1+X2+X3+Treatment,id=ID,data=M_data,corstr = "AR-M",Mv=1)
gee_model_2=gee(Outcome~X1+X2+X3+Treatment+Mediator,id=ID,data=Y_data,corstr = "AR-M",Mv=1)
gamma_gee=gee_model_2$coefficients[length(gee_model_2$coefficients)] 
effect_m_gee=gee_model_1$coefficients[length(gee_model_1$coefficients)]
indirect_effect_gee=gamma_gee*effect_m_gee
total_gee=indirect_effect_gee+gee_model_2$coefficients[length(gee_model_2$coefficients)-1] 
# Save results
save(pos_sample_size, direct_process,indirect_process,
     true_direct_effect_process,true_indirect_effect_process,
     true_total_effect_process,gamma_est,
     m_model_re,y_model_re,
     total_gee,effect_m_gee,gamma_gee,total_gee,
     total_effect,gamma,gamma_true,
     file=paste(sparsity.level,sample.size,i,"_simu_result.RData",sep="_"))
}