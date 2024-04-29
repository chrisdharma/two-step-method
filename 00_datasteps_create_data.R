load("cchsSNall.Rdata")

library(survey) 
library(tableone)
library(svrep)
library(dplyr)
library(tidyverse)

library(haven)
library(boot)
library(rsample)

#Some data manipulation
cchsSNall <- cchsSNall %>% mutate(employ = case_when(employ == 1 ~ 1,
                                                     employ %in% c(2,3,4) ~ 2,
                                                    employ %in% c(5,6) ~ 3))

names<-c("education","orient","immigration","existing_mh","consult_mh","alcohol","marijuana","smoke","insurance","illicit_drug","rural","poor_srmh",
         "have_pcp","employ","province","race","depress_binary","HIV_test12m","STI_test12m","lastcondom_use","survey","disclosure","loneliness_ucla_bin","mhsu_harmonized")
cchsSNall[,names]<-lapply(cchsSNall[,names],factor)

cchsSNall <- cchsSNall %>% 
  dplyr::select(c(existing_mh,depress_binary,loneliness_ucla_bin,
                  consult_mh,mhsu_harmonized,local_com_belong_bin,life_stress_bin,community_connection_lgb,
                  community_connection_msm,marijuana,poor_srmh,age_grp,income,education,orient,msm,employ,rural,race,
                  conversion_orientation,conversion_genderID,conversion_age,conversion_duration,conversion_num,outness,outness_age,
                  immigration,disclosure,survey,msm,year,ADM_RNO,WTS_M))

cchsSNall2 <- cchsSNall %>% drop_na(existing_mh,age_grp,income,education,orient,msm,employ,rural,race,immigration,disclosure,survey,msm,WTS_M)

###Now apply bootstrap, set up bootstrap data

setwd("/Users/christofferdharma/Documents/CCHS/AllDat")

bsw2015 <- read_sas("bsw2015.sas7bdat")
bsw2017 <- read_sas("bsw2017.sas7bdat")

#create bootstraps
B = 1000

cchs2015 <- cchsSNall2 %>% dplyr::filter(year == 2015)
cchs2017 <- cchsSNall2 %>% dplyr::filter(year == 2017)

sndat_all_bsw <- cchsSNall2 %>% dplyr::filter(survey == 1)

for (i in 1:B) {
  sndat_all_bsw[[as.symbol(paste0('BSW', i))]] <- 1
}

cchs2015_bsw <- merge(x = cchs2015, y = bsw2015, by = "ADM_RNO", all.x = TRUE)
cchs2017_bsw <- merge(x = cchs2017, y = bsw2017, by = "ADM_RNO", all.x = TRUE)

cchsbsw <- rbind(cchs2015_bsw,cchs2017_bsw)
cchsbsw <- cchsbsw  %>% dplyr::select(-c(fwgt))
nyear <- 2

divide <- function(x) (x/nyear)
cchsbsw <- cchsbsw %>% dplyr::mutate_at(.,vars(matches("BSW")), divide)
cchsbsw$WTS_M<-divide(cchsbsw$WTS_M)
allbsw <- rbind(cchsbsw,sndat_all_bsw)


