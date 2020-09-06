source("varcodebook.R")
#Removing missing covariates and mediation, 10216 remains out of 11525
Complete.data <-
  subset(gender.data,!apply(is.na(gender.data[, c(covariate.index, covariate.index.dsi,group_id,mediation.index[m], adverse.index[adv_id])]), 1, any))
#Truncate age at 18 
Complete.data<-subset(Complete.data,age_sample<=age.truncate.female)
#Order by sname and time
Complete.data<-Complete.data[order(Complete.data$sname,Complete.data$age_sample),]

#Exclude Those observed less than twice
Year.Obs.Info <-
  aggregate(
    Complete.data$age_sample,
    by = list(Complete.data$sname),
    FUN = function(x) {
      length(unique(x))
    }
  )
#Extract those have been observed at least three times
valid.id<-Year.Obs.Info$Group.1[Year.Obs.Info$x>=3]

length(valid.id)/length(unique(Year.Obs.Info$Group.1))
Complete.data<-subset(Complete.data,sname%in%valid.id)

#Log Transformation
Complete.data[,response.index]=log(Complete.data[,response.index])

#Transform Hydro Year into Categorical Variable
Complete.data$hydroyear=as.factor(Complete.data$hydroyear)
#Normalize the time to [0,1] Interval
range.y=range(Complete.data$age_sample)
Complete.data$t=(Complete.data$age_sample-range.y[1]+0.0001)/(range.y[2]-range.y[1]+0.0002)

Complete.data$treatment_factor=
  as.factor(Complete.data[,adverse.index[adv_id]])
