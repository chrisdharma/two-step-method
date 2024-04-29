library(survey) 
library(tableone)
library(svrep)
library(dplyr)
library(tidyverse)

library(haven)
library(boot)
library(rsample)

load("allbsw.Rdata")

M = 10
B = 100

allbsw <- allbsw %>% drop_na(age_grp,income,education,orient,msm,employ,rural,race,immigration,disclosure,survey,msm,WTS_M)

allmsm <- allbsw %>% filter(msm == 1)
allbsw$existing_mh <- as.numeric(as.character(allbsw$existing_mh))
allbsw$depress_binary <- as.numeric(as.character(allbsw$depress_binary))
allbsw$mhsu_harmonized <- as.numeric(as.character(allbsw$mhsu_harmonized))
allbsw$consult_mh <- as.numeric(as.character(allbsw$consult_mh))
allbsw$poor_srmh <- as.numeric(as.character(allbsw$poor_srmh))
allbsw$life_stress_bin <- as.numeric(as.character(allbsw$life_stress_bin))
allbsw$local_com_belong_bin <- as.numeric(as.character(allbsw$local_com_belong_bin))
allbsw$loneliness_ucla_bin <- as.numeric(as.character(allbsw$loneliness_ucla_bin))
allbsw$community_connection_msm <- as.numeric(as.character(allbsw$community_connection_msm))

allbsw$WTS_M <- as.numeric(as.character(allbsw$WTS_M))

######INITIAL STEP (1-3), which is the same for both approaches

#STEP 1. Calculate ALP on SN data
svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M, data = allmsm) ### 
Formula_fit.survey = as.formula("survey ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 

lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)    

beta.w = summary(lgtreg.w)$coeff[,1]   
allmsm$ps.w = lgtreg.w$fitted.values 

allmsm <- allmsm %>%
  mutate(ALP_weight = (1-ps.w) /ps.w,
         ALP_weight = ifelse(survey == 0,WTS_M,ALP_weight)) 

#Subset the data
SNall <- allmsm %>% dplyr::filter(survey == 1)
cchsall <- allbsw %>% dplyr::filter(survey == 0)

svyds.check.freq  = svydesign(ids=~1, weight = ~ WTS_M, data = cchsall) ### 
# svyCreateTableOne(vars=c("existing_mh"), data=svyds.check.freq)
# svyCreateTableOne(vars="existing_mh",factorVars="existing_mh",strata="msm",data=svyds.check.freq)

#STEP 2. Calculate pieces needed for Bayes formula
### We now want to predict probability of NOT disclosing, so create new var undisc
SNall$undisc <- relevel(as.factor(recode_factor(SNall$disclosure,'0'='1','1'='0')),ref="0")
svyds.ALP.SN = svydesign(ids=~1, weight = ~ ALP_weight, data = SNall)
Formula_fit.undisc = as.formula("undisc ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 
lgtreg.u = svyglm(Formula_fit.undisc, family= quasibinomial, design = svyds.ALP.SN)    

#Need two components, first, the odds of undiclosing, after applying ALP as shown here:
cchsall$odds<-exp(predict(lgtreg.u,newdata=cchsall,type="link")[1:nrow(cchsall)])

#Next, the probability of reporting positive in CCHS, conditional on X
svyds.cchs = svydesign(ids=~1, weight = ~ WTS_M, data = cchsall)
Formula_fit.msm = as.formula("msm ~ rural + age_grp + (age_grp**2) + income + education + race + employ + immigration + age_grp*employ") 
lgtreg.cchs.msm = svyglm(Formula_fit.msm, family= quasibinomial, design = svyds.cchs)    
cchsall$prep1<-predict(lgtreg.cchs.msm,newdata=cchsall,type="response")[1:nrow(cchsall)]

#Finally, need probability of reporting negative in CCHS, conditional on X, which is the inverse of the above
#So we can calculate pnewmsm as:
cchsall$pnewmsm <-cchsall$odds * ((cchsall$prep1) / (1-cchsall$prep1))

###3. Now apply MI on CCHS, function starts, to ensure we can run it on different outcomes
set.seed(365)

#Point estimate, working from CCHS side
point_est_mi_cchs <- function(inputdata,outcomevar) { 
msmimpute = vector(mode='list', length=M)
mhimp_cchs <- c()
newdat<-inputdata
newdat <- newdat %>% drop_na({{outcomevar}})
for (j in 1:M) {
  #If msm = 1, then 1, otherwise (0 or otherwise, it is a random number from the prob of pnewmsm)
  msmimpute[[j]]<-ifelse(newdat$msm==1,1,rbinom(size=1,n=nrow(newdat),prob=newdat$pnewmsm))
  newdat$msmimpute<-msmimpute[[j]]
  tempdat <- newdat %>% filter(msmimpute == 1)
  outcometemp <- tempdat %>% dplyr::select({{outcomevar}})
  weight <- tempdat %>% dplyr::select(c(WTS_M))
  mhimp_cchs[j]<- weighted.mean(outcometemp,weight)
  }
return(mean(mhimp_cchs))
# return(sum(newdat$WTS_M))
}
cchs_adj_existing_mh<-point_est_mi_cchs(inputdata=cchsall,outcomevar=existing_mh)
cchs_adj_mhsu<-point_est_mi_cchs(inputdata=cchsall,outcomevar=mhsu_harmonized)
cchs_adj_srmh<-point_est_mi_cchs(inputdata=cchsall,outcomevar=poor_srmh)
cchs_adj_stress<-point_est_mi_cchs(inputdata=cchsall,outcomevar=life_stress_bin)
cchs_adj_com<-point_est_mi_cchs(inputdata=cchsall,outcomevar=local_com_belong_bin)

####Now calculate the bootstrap input
bootstrap_cchs_matrix <- function(inputdata,outcomevar) {
  boot.estimate.cchs_adj <- matrix(nrow=B,ncol=M)
  newdat<-inputdata
  newdat <- newdat %>% drop_na({{outcomevar}})
for (j in 1:M) {
  newdat$msmimpute<-ifelse(newdat$msm==1,1,rbinom(size=1,n=nrow(newdat),prob=newdat$pnewmsm))
  tempdat <- newdat %>% filter(msmimpute == 1)
  imputed_estimate = vector()
  #bootstrap estimates to be selected randomly, every run, it will be different weights
  bootseq <- sample(seq(1:1000),B)
  for (i in 1:B) {
    bt_samples_cchs <- bootstraps(tempdat, times = B)
    x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
    variable = paste0("x.cchs$BSW",bootseq[i])
    outcometemp <- x.cchs %>% dplyr::select({{outcomevar}})
    weight <- eval(parse(text=variable))
    imputed_estimate[i]<- weighted.mean(outcometemp[[1]],weight)
  }
  boot.estimate.cchs_adj[,j] = imputed_estimate
  }
  return(boot.estimate.cchs_adj)
}
boot.estimate.cchs_adj_existing_mh<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=existing_mh)
boot.estimate.cchs_adj_depress_binary<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=depress_binary)
boot.estimate.cchs_adj_mhsu<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=mhsu_harmonized)
boot.estimate.cchs_adj_srmh<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=poor_srmh)
boot.estimate.cchs_adj_stress<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=life_stress_bin)
boot.estimate.cchs_adj_com<-bootstrap_cchs_matrix(inputdata=cchsall,outcomevar=local_com_belong_bin)

#Save the tables above for easier loading
# write.table(boot.estimate.cchs_adj_existing_mh, 'boot_mat_existing_mh.txt', col.names=NA)
# write.table(boot.estimate.cchs_adj_depress_binary, 'boot_mat_depress.txt', col.names=NA)
# write.table(boot.estimate.cchs_adj_mhsu, 'boot_mat_mhsu.txt', col.names=NA)
# write.table(boot.estimate.cchs_adj_srmh, 'boot_mat_srmh.txt', col.names=NA)
# write.table(boot.estimate.cchs_adj_stress, 'boot_mat_stress.txt', col.names=NA)
# write.table(boot.estimate.cchs_adj_com, 'boot_mat_com.txt', col.names=NA)

#Save them in a matrix

#Load the data in case you need them again
boot.estimate.cchs_adj_existing_mh <- read.delim('boot_mat_existing_mh.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.cchs_adj_mhsu <- read.delim('boot_mat_mhsu.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.cchs_adj_srmh <- read.delim('boot_mat_srmh.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.cchs_adj_stress <- read.delim('boot_mat_stress.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.cchs_adj_com <- read.delim('boot_mat_com.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)

bootstrap_ci <- function(boot.matrix,point.est) {
##Want to calculate M standard errors and M estimates from each column
var_impute_boot<-apply(boot.matrix,2,var)
mean_impute_boot<-apply(boot.matrix,2,mean)
#var is just the variance formula in R with n-1 (also B-1 for bootstrap, so this is correct) 
#margin = 2 for columns

#Then calculate standard MI CI var(MI)
var_impute <- sum(var_impute_boot)/M + (M+1)*sum((mean_impute_boot-point.est)^2)/(M*(M-1))
# And then take SE
sd_impute <- sqrt(var_impute)

#W is just the average variance from all imputations
W = mean(var_impute_boot)
V = sum((mean_impute_boot-point.est)^2) / (M-1)

#With t dist and R degrees of freedom
R <- (M-1)*(1 + ((M*W)/ ((M+1)*V)))^2
alpha = 0.05
t.score = qt(p=alpha/2, df=R,lower.tail=F)

lower = point.est - t.score * sd_impute
upper = point.est + t.score * sd_impute
return(c(lower=lower,est=point.est,upper=upper))
}

bootstrap_ci(boot.matrix=boot.estimate.cchs_adj_existing_mh,point.est=cchs_adj_existing_mh)
bootstrap_ci(boot.matrix=boot.estimate.cchs_adj_mhsu,point.est=cchs_adj_mhsu)
bootstrap_ci(boot.matrix=boot.estimate.cchs_adj_srmh,point.est=cchs_adj_srmh)
bootstrap_ci(boot.matrix=boot.estimate.cchs_adj_stress,point.est=cchs_adj_stress)
bootstrap_ci(boot.matrix=boot.estimate.cchs_adj_com,point.est=cchs_adj_com)

#####So now we're working from the Sex Now side
###Get the point estimate first, same thing, just need to create the function

point_est_mi_sn <- function(inputdatacchs,inputdatasn,outcomevar) { 
  cchsfunc<-inputdatacchs
  snfunc <- inputdatasn
  snfunc <- snfunc %>% drop_na({{outcomevar}})
  mhimp_SN <- c()
for (i in 1:M) {
  cchsfunc$newmsm<-ifelse(cchsfunc$msm==1,1,rbinom(size=1,n=nrow(cchsfunc),prob=cchsfunc$pnewmsm))
  cchsmsmnew <- cchsfunc %>% filter(newmsm == 1) %>% 
    dplyr::select(c({{outcomevar}},age_grp,income,education,orient,msm,employ,rural,race,immigration,disclosure,survey,msm,WTS_M))
  SNall2 <- snfunc %>% dplyr::select(c({{outcomevar}},age_grp,income,education,orient,msm,employ,rural,race,immigration,disclosure,survey,msm,WTS_M))
  newmsmdat <- rbind(cchsmsmnew,SNall2)
  
  svyds.wei  = svydesign(ids=~1, weight = ~ WTS_M, data = newmsmdat) ### 
  lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)    
  beta.w = summary(lgtreg.w)$coeff[,1]   
  newmsmdat$ps.w = lgtreg.w$fitted.values 
  
  newmsmdat <- newmsmdat %>%
    mutate(ALP_weight = (1-ps.w) /ps.w,
           ALP_weight = ifelse(survey == 0,WTS_M,ALP_weight))

  SNallnewALP <- newmsmdat %>% filter(survey == 1)
  outcometemp <- SNallnewALP %>% dplyr::select({{outcomevar}})
  outcome2 <- as.numeric(outcometemp[[1]])
  weight <- SNallnewALP %>% dplyr::select(c(ALP_weight))
  mhimp_SN[i]<-weighted.mean(x=outcome2,w=weight[[1]])
}
  return(mean(mhimp_SN))
  # return(nrow(snfunc))
}

sn_adj_existing_mh<-point_est_mi_sn(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=existing_mh) 
sn_adj_mhsu<-point_est_mi_sn(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=mhsu_harmonized) 
sn_adj_srmh<-point_est_mi_sn(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=poor_srmh) 
sn_adj_loneli<-point_est_mi_sn(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=loneliness_ucla_bin) 
sn_adj_msm_connect<-point_est_mi_sn(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=community_connection_msm) 

####Now calculate the bootstrap input, SN
bootstrap_sn_matrix <- function(inputdatacchs,inputdatasn,outcomevar) {
  boot.estimate.sn_adj <- matrix(nrow=B,ncol=M)
  cchsfunc<-inputdatacchs
  snfunc <- inputdatasn
  snfunc <- snfunc %>%
    dplyr::select(-c(ps.w,ALP_weight,undisc)) %>% 
    drop_na({{outcomevar}})
for (j in 1:M) {
  cchsfunc$newmsm<-ifelse(cchsfunc$msm==1,1,rbinom(size=1,n=nrow(cchsfunc),prob=cchsfunc$pnewmsm))
  cchsmsmnew <- cchsfunc %>% filter(newmsm == 1) %>% 
    dplyr::select(-c(newmsm,pnewmsm,odds,prep1))
  #bootstrap estimates to be selected randomly, every run, it will be different weights
  bootseq <- sample(seq(1:1000),B)
  bt_samples_cchs <- bootstraps(cchsmsmnew, times = B)
  bt_samples_sn <- bootstraps(snfunc, times = B)
  imputed_estimate<-vector()
  for (i in 1:B) {
    x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
    x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
    x.all<-rbind(x.cchs,x.sn)
    variable = paste0("BSW",bootseq[i])
    svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = x.all) ### 
    lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)   
    beta.w = summary(lgtreg.w)$coeff[,1]   
    x.all$ps.w = lgtreg.w$fitted.values 
    SN.all.ALP.last <- x.all %>% dplyr::filter(survey == 1)
    SN.all.ALP.last <- SN.all.ALP.last %>%
      mutate(ALP_weight = (1-ps.w) /ps.w)
    outcometemp <- SN.all.ALP.last %>% dplyr::select({{outcomevar}})
    outcome2 <- as.numeric(outcometemp[[1]])
    weight <- SN.all.ALP.last %>% dplyr::select(c(ALP_weight))
    imputed_estimate[i]<-weighted.mean(x=outcome2,w=weight[[1]])
    }
  boot.estimate.sn_adj[,j] = imputed_estimate
  }
return(boot.estimate.sn_adj)  
}

boot.estimate.sn_adj_existing_mh<-bootstrap_sn_matrix(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=existing_mh)
boot.estimate.sn_adj_mhsu<-bootstrap_sn_matrix(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=mhsu_harmonized)
boot.estimate.sn_adj_srmh<-bootstrap_sn_matrix(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=poor_srmh)
boot.estimate.sn_adj_loneli<-bootstrap_sn_matrix(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=loneliness_ucla_bin)
boot.estimate.sn_adj_msm_connect<-bootstrap_sn_matrix(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=community_connection_msm)

setwd("/Users/christofferdharma/Documents/ThesisDat/boot_matrix/sn")

write.table(boot.estimate.sn_adj_existing_mh, 'boot_mat_existing_mh.txt', col.names=NA)
write.table(boot.estimate.sn_adj_mhsu, 'boot_mat_mhsu.txt', col.names=NA)
write.table(boot.estimate.sn_adj_srmh, 'boot_mat_srmh.txt', col.names=NA)
write.table(boot.estimate.sn_adj_loneli, 'boot_mat_loneli.txt', col.names=NA)
write.table(boot.estimate.sn_adj_msm_connect, 'boot_mat_msmcom.txt', col.names=NA)

#Load the data in case you need them again
boot.estimate.sn_adj_existing_mh <- read.delim('boot_mat_existing_mh.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.sn_adj_mhsu <- read.delim('boot_mat_mhsu.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.sn_adj_srmh <- read.delim('boot_mat_srmh.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.sn_adj_loneli <- read.delim('boot_mat_loneli.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)
boot.estimate.sn_adj_msm_connect <- read.delim('boot_mat_msmcom.txt', header = TRUE, sep = ' ', dec = ".",row.names=1)

bootstrap_ci(boot.matrix=boot.estimate.sn_adj_existing_mh,point.est=sn_adj_existing_mh)
bootstrap_ci(boot.matrix=boot.estimate.sn_adj_mhsu,point.est=sn_adj_mhsu)
bootstrap_ci(boot.matrix=boot.estimate.sn_adj_srmh,point.est=sn_adj_srmh)
bootstrap_ci(boot.matrix=boot.estimate.sn_adj_loneli,point.est=sn_adj_loneli)
bootstrap_ci(boot.matrix=boot.estimate.sn_adj_msm_connect,point.est=sn_adj_msm_connect)

###Now calculate unadjusted cis
cchsall_msm <- cchsall %>% dplyr::filter(msm == 1)

boot.cchs_unadj <- function(outcomevar) {
  newdat <- cchsall_msm %>% drop_na({{outcomevar}})
  outcomeall <- newdat %>% dplyr::select({{outcomevar}})
  estimate<- weighted.mean(outcomeall[[1]],newdat$WTS_M)
  bootseq <- sample(seq(1:1000),B)
  bs_estimates<-c()
  for (i in 1:B) {
  bt_samples_cchs <- bootstraps(newdat, times = B)
  x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
  variable = paste0("x.cchs$BSW",bootseq[i])
  outcometemp <- x.cchs %>% dplyr::select({{outcomevar}})
  weight <- eval(parse(text=variable))
  bs_estimates[i]<- weighted.mean(outcometemp[[1]],weight)
  }
  return(c(lower=quantile(bs_estimates,.025),est=estimate,upper=quantile(bs_estimates,.975)))
}
boot.cchs_unadj(outcomevar=existing_mh)
boot.cchs_unadj(outcomevar=mhsu_harmonized)
boot.cchs_unadj(outcomevar=poor_srmh)
boot.cchs_unadj(outcomevar=life_stress_bin)
boot.cchs_unadj(outcomevar=local_com_belong_bin)

boot.cchs_unweighted <- function(outcomevar) {
  newdat <- cchsall_msm %>% drop_na({{outcomevar}})
  outcomeall <- newdat %>% dplyr::select({{outcomevar}})
  estimate<- mean(outcomeall[[1]])
  bootseq <- sample(seq(1:1000),B)
  bs_estimates<-c()
  for (i in 1:B) {
    bt_samples_cchs <- bootstraps(newdat, times = B)
    x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
    outcometemp <- x.cchs %>% dplyr::select({{outcomevar}})
    bs_estimates[i]<- mean(outcometemp[[1]])
  }
  return(c(lower=quantile(bs_estimates,.025),est=estimate,upper=quantile(bs_estimates,.975)))
}
boot.cchs_unweighted(outcomevar=existing_mh)
boot.cchs_unweighted(outcomevar=mhsu_harmonized)
boot.cchs_unweighted(outcomevar=poor_srmh)
boot.cchs_unweighted(outcomevar=life_stress_bin)
boot.cchs_unweighted(outcomevar=local_com_belong_bin)

boot.sn_unadj <- function(outcomevar) {
  newdat <- SNall %>% drop_na({{outcomevar}})
  outcomeall <- newdat %>% dplyr::select({{outcomevar}})
  estimate<- prop.table(table(outcomeall[[1]]))[[2]]
  bs_estimates<-c()
  for (i in 1:B) {
    bt_samples_sn <- bootstraps(newdat, times = B)
    x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
    outcometemp <- x.sn %>% dplyr::select({{outcomevar}})
    bs_estimates[i]<-prop.table(table(outcometemp[[1]]))[[2]]
  }
  return(c(lower=quantile(bs_estimates,.025),est=estimate,upper=quantile(bs_estimates,.975)))
}
boot.sn_unadj(outcomevar=existing_mh)
boot.sn_unadj(outcomevar=mhsu_harmonized)
boot.sn_unadj(outcomevar=poor_srmh)
boot.sn_unadj(outcomevar=loneliness_ucla_bin)
boot.sn_unadj(outcomevar=community_connection_msm)

##ALP estimates for bootstrap, no MI

bootstrap_sn_ALP <- function(inputdatacchs,inputdatasn,outcomevar) {
  #To calculate point estimates, we already have ALP_weight, so just use them
  sntemp <- inputdatasn %>% drop_na({{outcomevar}})
  outcometemp <- sntemp %>% dplyr::select({{outcomevar}}) 
  outcome2 <- as.numeric(as.character(outcometemp[[1]]))
  weight <- sntemp %>% dplyr::select(c(ALP_weight))
  estimate<-weighted.mean(x=outcome2,w=weight[[1]])
  
  snfunc <- inputdatasn %>%
    dplyr::select(-c(ps.w,ALP_weight,undisc)) %>% 
    drop_na({{outcomevar}})
  
  #The CCHS people here needs to be all SMM
  cchsfunc <- inputdatacchs %>%
    dplyr::select(-c(odds,prep1,pnewmsm)) %>% 
    filter(msm == 1)
  
  bt_samples_cchs <- bootstraps(cchsfunc, times = B)
  bt_samples_sn <- bootstraps(snfunc, times = B)
  bootseq <- sample(seq(1:1000),B)
  imputed_estimate<-vector()

  for (i in 1:B) {
    
  x.cchs<-analysis(bt_samples_cchs$splits[[i]]) %>% as_tibble()
  x.sn<-analysis(bt_samples_sn$splits[[i]]) %>% as_tibble()
  x.all<-rbind(x.cchs,x.sn)
  variable = paste0("x.all$BSW",bootseq[i])
  svyds.wei  = svydesign(ids=~1, weight = ~ eval(parse(text=variable)), data = x.all) ### 
  lgtreg.w = svyglm(Formula_fit.survey, family= quasibinomial, design = svyds.wei)   
  beta.w = summary(lgtreg.w)$coeff[,1]   
  x.all$ps.w = lgtreg.w$fitted.values 
  SNonly <- x.all %>% dplyr::filter(survey == 1)
  SNonly <- SNonly %>%
    mutate(ALP_weight = (1-ps.w) /ps.w)
  sntemp <- SNonly %>% drop_na({{outcomevar}})
  outcometemp <- sntemp %>% dplyr::select({{outcomevar}}) 
  outcome2 <- as.numeric(as.character(outcometemp[[1]]))
  outcome2 <- as.numeric(as.character(outcometemp[[1]]))
  weight <- sntemp %>% dplyr::select(c(ALP_weight))
  imputed_estimate[i]<-weighted.mean(x=outcome2,w=weight[[1]])
  }
  return(c(lower=quantile(imputed_estimate,.025),est=estimate,upper=quantile(imputed_estimate,.975)))
}

bootstrap_sn_ALP(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=existing_mh)
bootstrap_sn_ALP(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=mhsu_harmonized)
bootstrap_sn_ALP(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=poor_srmh)
bootstrap_sn_ALP(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=loneliness_ucla_bin)
bootstrap_sn_ALP(inputdatacchs=cchsall,inputdatasn=SNall,outcomevar=community_connection_msm)

