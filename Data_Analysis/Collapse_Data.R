###Load Data
#rm(list=ls())
adverse.index <-
  c(
    "adv_rain",
    "adv_sib",
    "adv_density",
    "adv_mom",
    "adv_mom_rank",
    "adv_mom_dsi",
    "adv_cumulative",
    "adv_01",
    "adv_12",
    "adv_02"
  )

adverse.name <-
  c(
    "Rain Adversity",
    "Sibling Adversity",
    "Density Adversity",
    "Maternal Death",
    "Maternal Rank",
    "Maternal DSI",
    "Cumulative Adversity",
    "Cumulative 0/1 Compared",
    "Cumulative 1/2 Compared",
    "Cumulative 0/2 Compared"
  )


covariate.index.female = c(
  "state",
  "n_adults",
  "n_adults_squared",
  "mean_tmax_30d",
  "season",
  "rain_anom_3mo",
  #  "hyb_score",
  #  "hydroyear",
  "n_adult_mat_kin"
)

covariate.index.male = c(
  "n_adults",
  "n_adults_squared",
  "mean_tmax_30d",
  "season",
  "rain_anom_3mo",
  #  "hyb_score",
  #"hydroyear",
  "n_adult_mat_kin"
)

mediation.index <- c("DSI_F", "DSI_M", "SCI_F", "SCI_M", "proprank")

#setwd("~/Downloads/Second Year/FPCA/wrap/newdata")
load("ELADSIGC data for Shuxi 9.22.2019.RData")
#Divide by Gender
if(f==1){
  gender.data=subset(merged_adv_gcdata,sex=='F')
}else{
  gender.data=subset(merged_adv_gcdata,sex=='M')
}


#Create Additional Treatment Indicator
#0/1 Compared
gender.data$adv_01 = FALSE
#Additional Condition
#gender.data$adv_01[gender.data$adv_cumulative > 0&gender.data$adv_mom_rank==0] = TRUE
#gender.data$adv_01[gender.data$adv_cumulative > 0&gender.data$adv_mom_rank==0] = TRUE
gender.data$adv_01[gender.data$adv_cumulative > 0&gender.data$adv_rain==0&gender.data$adv_mom_rank==0] = TRUE


#1/2 Compared
gender.data$adv_12 = NA
gender.data$adv_12[gender.data$adv_cumulative == 1&gender.data$adv_rain==0&gender.data$adv_mom_rank==0] = FALSE
gender.data$adv_12[gender.data$adv_cumulative > 1&gender.data$adv_rain==0&gender.data$adv_mom_rank==0] = TRUE

#1/2 Compared
gender.data$adv_02 = NA
gender.data$adv_02[gender.data$adv_cumulative == 0] = FALSE
gender.data$adv_02[gender.data$adv_cumulative > 2] = TRUE


###Collapse
cov_continu_data=aggregate(gender.data[,c("gc",mediation.index,c("adv_cumulative","n_adults","n_adults_squared",
                                                                 "mean_tmax_30d","rain_anom_3mo",  
                                                                 "n_adult_mat_kin","avg_mat_kin",
                                                                 "percent_days_with_infant",'grp',
                                                                 "percent_days_cycling","hyb_score",
                                                                 "avg_density","avg_age","hydroyear"))],
                           by=list(gender.data$age_sample,gender.data$sname),
                           mean)
adv_data = aggregate(gender.data[, adverse.index[-7]],
                     by = list(gender.data$age_sample, gender.data$sname),
                     any)
cat_data=aggregate(gender.data[,c('season','state')],by = list(gender.data$age_sample,gender.data$sname),
                   FUN=function(x){sample(x,1)})
gender.data=cbind(cov_continu_data,adv_data[,-c(1,2)],cat_data[,-c(1,2)])

colnames(gender.data)[c(1,2)]=c("age_sample","sname")
