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
  "proprank"
# "hyb_score"
#  "hydroyear",
#  "n_adult_mat_kin"
)

covariate.index.male = c(
  "n_adults",
  "n_adults_squared",
  "mean_tmax_30d",
  "season",
  "rain_anom_3mo"
  #"proprank",
  #  "hyb_score"
  #"hydroyear",
  #n_adult_mat_kin"
)

#Determine the covariates for the mediators
if(f==1){
  if(m==1||m==3){
covariate.index.dsi=c(
  #"percent_days_cycling",
  "percent_days_with_infant",
  "avg_mat_kin",
  "avg_density",
  "avg_age"
#  "hydroyear"
)
}else if(m==2||m==4)
{
  covariate.index.dsi=c(
    "percent_days_cycling",
    # "percent_days_with_infant",
    "avg_density",
    "avg_age"
    #  "hydroyear"
  )
}else{  covariate.index.dsi=c(
"n_adults"
  #  "hydroyear"
)}}else{
  if(m==1||m==3){
    covariate.index.dsi=c(
      "avg_mat_kin",
      "avg_density",
      "avg_age"
      #  "hydroyear"
    )
  }else{
    covariate.index.dsi=c(
      "n_adults"
      #  "hydroyear"
    )
    covariate.index.female = c(
      "state",
      "n_adults",
      "n_adults_squared",
      "mean_tmax_30d",
      "season",
      "rain_anom_3mo"
      #"proprank",
      #"hyb_score"
      #  "hydroyear",
#      "n_adult_mat_kin"
    )
    
  }
  
}


mediation.index <- c("DSI_F", "DSI_M", "SCI_F", "SCI_M", "proprank")
response.index <- c("gc")
group_id=c("grp")

if(f==1)
{
  covariate.index=covariate.index.female
}else{
  covariate.index=covariate.index.male
}
covariate.collect.index=union(covariate.index.dsi,covariate.index)



