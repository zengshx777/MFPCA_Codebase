##Main Script to Run to analysis 
rm(list=ls())
args=commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#setwd("~/Downloads/Second Year/FPCA/wrap/newdata")


#mediation index
# #m=as.numeric(args[1])

# #Adversity index
# #adv_id=as.numeric(args[2])

# #sex index
# #f=as.numeric(args[3])

##Grids for the orthogonality condition
Grid.num=40
work_grid=seq(0,1,length=Grid.num)

# m=1
# adv_id=10
# f=1
# age.truncate.female=17
# l_m=11;K_m=7
# l_y=11;K_y=7
##MCMC Time
simu_time=20000;

result_sum=NULL
source("Collapse_Data.R")
# for (adv_id in 1:9)
# {
tryCatch({
source("Clean_Data.R")
age_grids=work_grid*(range.y[2]-range.y[1]+0.0002)+range.y[1]-0.0001
##Prior Parameter, if we assume the causal parameters have a MGP
tau_mgp=1
##If we assume individual precision paramter on noise term
kappa_sample=0
##If have individual precision parameter on principal score
phi_sample=1
##If enforcing deterministic increasing
truncate=0

##Whether include time main effect 1 for no/0 for yes
simple=1;


l=l_m;#Number of Basis constructing the PC, larger than K
K=K_m#Number of PC
bandwidth=4 #Bandwidth for exponential kernel
#Adaptive Shrinkage on PC
shrink=1
source("M_stage.R")

l=l_y;K=K_y#Number of PC
bandwidth=4
###Sampler for the outcome
#source("YStage_New.R")
source("Y_stage.R")
save(Y_MCMC_Result,
     M_MCMC_Result,
     Complete.data,
     l_y,l_m,K_y,K_m,age_grids)
# 
# #pdf(file=paste(adverse.index[adv_id],mediation.index[m],f,"Decom.pdf",sep="_"),height=8,width=12)
# plot(total_effect_mean,ylim=range(total_effect),type='l',col="blue",
#      main=paste("Effect Decomposition\n",adverse.name[adv_id]),ylab="Effect on GC")
# lines(total_effect_down,lty=2)
# lines(total_effect_up,lty=2)
# abline(h=0)
# 
# lines(indirect_process_mean,col="red")
# lines(indirect_process_up,lty=2)
# lines(indirect_process_down,lty=2)
# abline(h=0)
# 
# legend("topright",legend=c("Total Effect","Mediation Effect"),col=c("blue","red"),
#        lty=1)
#dev.off()

result_sum=rbind(result_sum,c(mean(mediator_effect_mean),
             mean(mediator_effect_down),
             mean(mediator_effect_up),
             mean(gamma),
             quantile(gamma,0.025),
             quantile(gamma,0.975),
             mean(indirect_process_mean),
             mean(indirect_process_down),
             mean(indirect_process_up),
             mean(total_effect_mean),
             mean(total_effect_down),
             mean(total_effect_up)))
},error=function(e){
  print(e)
  result_sum=rbind(result_sum,rep(NA,12))
})
print(adverse.name[adv_id])
#}
setwd("Results/")
save(result_sum,mediator_effect_mean,mediator_effect_down,mediator_effect_up,
     gamma,indirect_process_mean,indirect_process_down,indirect_process_up,
     total_effect_mean,total_effect_down,total_effect_up,
     file=paste(mediation.index[m],adv_id,f,"Collect_withHybScore.RData",sep="_"))
