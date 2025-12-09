#### Begin - Set Working Directory ####
## quick script for cleaning environmental data
setwd("C:/Users/trevo/Dropbox/My PC (LAPTOP-GI7LHD15)/Documents/GitHub/Fire-Retartand-Veg-Soils") ## set working directory to the data folder


trtdata <- read.csv("./data/plotMetadata.csv")
sc_env<- read.csv("./data/StoneCanyon_ENV.csv")
quarry_env<- read.csv("./data/Quarry_ENV.csv")

colnames(trtdata)
env <- data.frame(plot = trtdata$plot,
                  site = trtdata$site,
                  sev = trtdata$sev,
                  trt = trtdata$trt,
                  elev = NA,
                  slope = NA,
                  tpi = NA,
                  aspect = NA)
env$elev[match(sc_env$Plot,env$plot)] <- sc_env$Elevation
env$slope[match(sc_env$Plot,env$plot)] <- sc_env$Slope
env$tpi[match(sc_env$Plot,env$plot)] <- sc_env$TPI
env$aspect[match(sc_env$Plot,env$plot)] <- sc_env$aspect

env$elev[match(quarry_env$GPX,env$plot)] <- quarry_env$Elevation
env$slope[match(quarry_env$GPX,env$plot)] <- quarry_env$Slope
env$tpi[match(quarry_env$GPX,env$plot)] <- quarry_env$TPI
env$aspect[match(quarry_env$GPX,env$plot)] <- quarry_env$aspect

write.csv(env, "./data/ENV_Data.csv")
rm(list = ls()) ## cleaning global env
gc()
