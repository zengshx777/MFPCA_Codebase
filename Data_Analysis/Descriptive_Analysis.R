###Descriptive Analysis Part

rm(list=ls())
# setwd("C:/Users/Shuxi ZENG/Dropbox/Third Year/FPCA_New/FPCA_1005")
load("ELADSIGC data for Shuxi 9.22.2019.RData")

f = 1
age.truncate.female = 18
source("Collapse_Data.R")

covariate_descriptive = c(
  "age_sample",
  "n_adults",
  "avg_mat_kin",
  "percent_days_cycling",
  "percent_days_with_infant",
  "mean_tmax_30d",
  "rain_anom_3mo",
  "proprank",
  "hyb_score"
)
table_num = NULL

#Number of individuals
for (m in 1:5)
{
  for (adv_id in (1:9))
  {
    source("Clean_Data.R")
    
    
    table_num = rbind(table_num, c(
      m,
      adv_id,
      length(unique(Complete.data$sname)),
      length(unique(Complete.data$grp)),
      length(Complete.data$sname)
    ))
  }
}
colnames(table_num) = c(
  "Mediator Index",
  "Adversity Index",
  "Number of Individuals",
  "Number of Social Groups",
  "Sample Size"
)

write.csv(table_num, file = "number_summary.csv")


f=1;m=1;adv_id=1
source("Clean_Data_0923.R")
summary_cov=t(apply(
  Complete.data[, covariate_descriptive],
  2,
  FUN = function(x) {
    c(mean(x[!is.na(x)]), median(x[!is.na(x)]), sd(x[!is.na(x)]),range(x[!is.na(x)]))
  }
))
colnames(summary_cov)=c("mean","median","sd","min","max")
write.csv(summary_cov,file="covariate_descriptive.csv")

write.csv(c(table(Complete.data$state),
table(Complete.data$hydroyear),
table(Complete.data$season)),file="categorical_variable_descriptive.csv")


##Adversity Counts Histogram
hist(
  gender.data$adv_cumulative,
  breaks = rep(0:5, each = 2) + c(-0.4, 0.4),
  main = "Histogram of Cumulative Adversity Counts",
  freq = TRUE,
  xlab = "Cumulative Adversities",
  ylab = "Sample Frequency"
)

hist(aggregate(Complete.data$adv_cumulative,by=list(gender.data$sname),mean)$x,
     breaks = rep(0:5, each = 2) + c(-0.4, 0.4),
     main = "Histogram of Cumulative Adversity Counts",
     freq = TRUE,
     xlab = "Cumulative Adversities",
     ylab = "Count of Females")

write.csv(t(unlist(apply(
  Complete.data[, adverse.index[-7]],
  2,
  FUN = function(x) {
    table(x)
  }
))),file="adv_sample_counts.csv")

individual_accum=aggregate(Complete.data$adv_cumulative,by=list(Complete.data$sname),mean)$x
adv_individual_summary=c(mean(individual_accum),median(individual_accum),
  sd(individual_accum),range(individual_accum))
names(adv_individual_summary)=c("mean","median","sd","min","max")
write.csv(adv_individual_summary,file="accumulative_individual_summary.csv")

write.csv(t(unlist(apply(
  aggregate(
    Complete.data[, adverse.index[-7]],
    by = list(Complete.data$sname),
    FUN = function(x) {
      mean(x)
    }
  )[, -1], 2, table
))),file="adv_individual_counts.csv")
#Number of GC samples
gc_obs_counts=aggregate(Complete.data$gc,by=list(Complete.data$sname),length)$x
hist(gc_obs_counts,
     breaks = rep(c(0,50,100,150,200,250,300), each = 2)+c(-25,25),
     freq=TRUE,
     xlab="Number of GC Observations",
     ylab="Sample Counts",
     main="Histogram of GC Observations")

gc_counts=c(mean(gc_obs_counts),median(gc_obs_counts),sd(gc_obs_counts),range(gc_obs_counts))
names(gc_counts)=c("mean","median","sd","min","max")
write.csv(gc_counts,file="gc_sample_counts.csv")


###Mediator Description
mediator_counts_table=NULL
for (m in 1:5)
{
  source("Clean_Data_0923.R")
  unique_counts=aggregate(Complete.data$hydroyear,by=list(Complete.data$sname),FUN=function(x){length(unique(x))})$x
  mediator_counts_table=rbind(mediator_counts_table,
        c(mean(unique_counts),median(unique_counts),sd(unique_counts),range(unique_counts),sum(unique_counts) ))
}
row.names(mediator_counts_table)=mediation.index
colnames(mediator_counts_table)=c("mean","median","sd","min","max","total")
write.csv(mediator_counts_table,file="mediator_obs_counts.csv")