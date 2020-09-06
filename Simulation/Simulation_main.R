#Illustration
source("Simulation_Function.R")
#Assumption 1,2,3 Valid
Simu_Wrap=Simulation_Data(nsample_m_mean=15,
                          nsample_y_mean=15,
                          nsample_size=200,
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
save(pos_sample_size, direct_process,indirect_process,
     total_effect,gamma,file="S1.RData")

true_direct_effect_process = mean(y_1_1-y_0_1)
true_indirect_effect_process = mean(y_0_1-y_0_0)
true_total_effect_process = mean(y_1_1-y_0_0)
gamma_true =1 
direct_effect_process=apply(simplify2array(direct_process),1:2,mean)
indirect_effect_process=apply(simplify2array(indirect_process),1:2,mean)
total_effect_process=apply(simplify2array(total_effect),1:2,mean)
gamma_est = mean(gamma)

(mean(indirect_effect_process)-mean(true_indirect_effect_process))

## RANDOM
library(mgcv)
model_1=gam(Mediator~X1+X2+X3+Treatment+s(Time)+s(ID,bs="re") ,data=M_data)
Y_data$Mediator = M_data$Mediator
model_2=gam(Outcome~X1+X2+X3+Mediator+Treatment+s(Time)+s(ID,bs="re") ,data=Y_data)



## GEE
library(gee)
gee_model_1=gee(Mediator~X1+X2+X3+Treatment,id=ID,data=M_data,corstr = "AR-M",Mv=1)
gee_model_2=gee(Outcome~X1+X2+X3+Treatment+Mediator,id=ID,data=Y_data,corstr = "AR-M",Mv=1)
gamma_gee=gee_model_2$coefficients[length(gee_model_2$coefficients)] 
effect_m_gee=gee_model_1$coefficients[length(gee_model_1$coefficients)]
indirect_effect_gee=gamma_gee*effect_m_gee
total_gee=indirect_effect_gee+gee_model_2$coefficients[length(gee_model_2$coefficients)-1] 

