#### Begin - Loading Dependencies, Set Working Directory ####
library(dplyr) ## useful for transforming data
library(stringr) ## needed to find patterns in character strings
library(lme4) ## glmm package
library(vegan) ## multidimentional vegetation analysis package
setwd("C:/Users/trevo/Dropbox/My PC (LAPTOP-GI7LHD15)/Documents/GitHub/Fire-Retartand-Veg-Soils") ## set working directory to the main GitHub folder

#### Reading in Data ####
cover <- read.csv("./data/FieldData(comm).csv") ## field data on plant cover
rich <- read.csv("./data/FieldData(rich).csv") ## field data on species richness
env <- read.csv("./data/ENV_Data.csv") ## environmental data from script 1-Clean Environmental Data.R
sp.info <- read.csv("./data/IntroducedStatusKey.csv") ## species information on introduced/native, functional group, etc. 
soils <- read.csv("./data/CarterFR_analysis.csv")

#### cleaning data ####
## matching the plot names, a bit of inconsistency
cover$plot <- gsub("_", "", cover$plot)
rich$plot <- gsub("_", "", rich$plot)
env$plot <- gsub("_", "", env$plot)
soils$Field.ID <- gsub("_", "", soils$Field.ID)

soils <- soils[complete.cases(soils),] ## soils is relatively clean already
soils <- soils[order(soils$Location,soils$Treatment,soils$Burn.Severity),] ## ordering the data 

`%notin%` <- purrr::negate(`%in%`) ## creating a function to negate %in% , makes cleaning easier
env$plot[env$plot[order(env$plot)] %notin% soils$Field.ID[order(soils$Field.ID)]]
## this looks correct, MFR05 and SCCON11 were not visited for sampling

env <- env[match(soils$Field.ID, env$plot),] ## correct number of plots

length(unique(sp.info$code)) ## information that is most clean (no duplicates); 231 species
length(unique(cover$code)) ## 166 species
length(unique(rich$code)) ## 253 species

unique(cover$code[cover$code %notin% sp.info$code]) ## finding observations not in the sp.info df
cover <- cover[cover$code != "ROCK" & 
                 cover$code != "LITTER" & 
                 cover$code != "LITTER " &
                 cover$code != "LITTER`" &
                 cover$code != "SOIL",] ## removing non-species codes
## manually fixing spelling errors
cover$code[cover$code == "ANDGER "] <- "ANDGER"
cover$code[cover$code == "CAREX"] <- "CAREX SP."
cover$code[cover$code == "ESCVIR"] <- "ECHVIR" 
cover$code[cover$code == "UNK 04" | cover$code == "UNK 05" | cover$code == "UNK 14"] <- "UNABLETOID"
cover$code[cover$code == "BOUSTR"] <- "BOESTR" 
cover$code[cover$code == "LATALN"] <- "LATLAN" 
cover$code[cover$code == "POLDEL"] <- "POLDOU" 
cover$code[cover$code == "CAREX "] <- "CAREX SP."
cover$code[cover$code == "AQUCER"] <- "AQUCOE"
cover$code[cover$code == "MEHRAN"] <- "MAHREP"

## not sure what these codes were supposed to be so I am removing them
cover <- cover[cover$code != "HIEALB",]  # no clue, removing for now (will look at data on campus)
cover <- cover[cover$code != "LONUTA",]  # no clue, removing
unique(cover$code[cover$code %notin% sp.info$code]) ## now the only missing code is the UNABLETOID code, which makes sense
length(unique(cover$code)) ## 150 species in the cover data

## repeating the process for the species richness df
unique(rich$code[rich$code %notin% sp.info$code]) ## mostly unknowns as errors
rich$code[rich$code == "UNK10" | 
            rich$code == "UNK11" | 
            rich$code == "UNK 11" | 
            rich$code == "UNK 09" | 
            rich$code == "UNK 04" | 
            rich$code == "UNK 02" | 
            rich$code == "UNK6" | 
            rich$code == "UNK13" | 
            rich$code == "UNK12" | 
            rich$code == "UNK14" | 
            rich$code == "UNK15" | 
            rich$code == "UNK A"] <- "UNABLETOID"
rich$code[rich$code == "CAREX sp."] <- "CAREX SP."
rich$code[rich$code == "HUEPAR"] <- "HEUPAR"
rich$code[rich$code == "ARTDRA "] <- "ARTDRA"
rich$code[rich$code == "ESCVIR"] <- "ECHVIR"
rich$code[rich$code == "CIRSIUM SP."] ## not doing anything with the cirsium genus at the moment
rich$code[rich$code == "RHUARO"] <- "RHUTRI"
rich$code[rich$code == "OPUPOL "] <- "OPUPOL"
rich$code[rich$code == "BOUSTR"] <- "BOESTR"
rich$code[rich$code == "HETHIR"] <- "HETVIL"
rich <- rich[rich$code != "GILCUD",] ## no clue, removing for now

unique(rich$code[rich$code %notin% sp.info$code])
length(unique(rich$code)) ## 233 species

#### Average Species Cover per Plot ####
table(cover$plot) ## looking at plots and double checkking
length(unique(cover$plot)) ## No plots missing
table(rich$plot)
length(unique(rich$plot)) ## No plots missing

length(unique(cover$code))
str(cover)
cover$cov <- as.numeric(cover$cov)
cover[is.na(cover$cov),] ## will have to fix this later

trans.avg <- function(x){
  transavg <- (sum(x))/3
  return(transavg)
} ## Custom function to take the average on each transect. 
## Sum is divided by 3 because there were 3 subplots and 0 values were not recorded.
plot.avg <- function(x){
  plotavg <- (sum(x))/3
  return(plotavg)
} ## Custom function to take the average on each plot
## Sum is divided by 3 because there were 3 transects and 0 values were not recorded.

cover.sum <- cover %>% 
  group_by(plot, transect, code) %>% 
  summarise(transavg = trans.avg(cov)) %>%
  group_by(plot, code) %>%
  summarise(plotavg = plot.avg(transavg)) ## summarizing by transect then by plot

str(cover.sum)
cover.sum <- as.data.frame(cover.sum)
length(unique(cover.sum$plot)) ## 108 plots (this is correct)

#### Creating Community Data (Cover Estimates + Rare Species only found as richness) ####
rich
rich <- rbind(rich[,c(1,2)], cover[,c(1,4)]) ## adding the rich and cover data together, just to make sure we didn't miss any species in any plot
rich <- rich[order(rich$plot, rich$code),] ## ordering the dataframe
rich <- rich[!duplicated(rich),] ## removing duplicates (same species and same plot, but not same species different plot)
table(rich$code)
rich$plotavg <- 0

comm <- rbind(cover.sum,rich) ## creating a community dataframe that has both cover and richness
comm <- comm[order(comm$plot, comm$code),] ## ordering the community data

comm.sum <- NA
vec <- unique(comm$plot)
for(i in 1:length(vec)){
  tmp <- comm[comm$plot == vec[i],]
  tmp <- tmp[!duplicated(tmp$code),]
  comm.sum <- rbind(comm.sum,tmp)
  print(i/length(vec))
} ## removing the richness (with cover of 0) for species that were present in the cover estimates
comm.sum <- comm.sum[-1,] ## removing the first row (error)
comm.sum <- comm.sum[order(comm.sum$plot, comm.sum$code,comm.sum$plotavg),] ## ordering the data

comm.sum.long <- comm.sum ## making a second community dataframe that I can transform and use for further analysis
comm.sum <- reshape(comm.sum, idvar = "plot", timevar = "code", direction = "wide") ## reshaping the data
rownames(comm.sum) <- comm.sum$plot ## renaming rows to be PlotID
comm.sum$plot <- NULL ## removing non species column
colnames(comm.sum) <- sub("plotavg.", "", colnames(comm.sum)) ## changing column names to be species codes

PlotCov <- as.data.frame(apply(comm.sum, 1, sum, na.rm = TRUE)) ## getting the total cover per plot
colnames(PlotCov) <- "TotalPlotCover" ## renaming the column for later

SpeciesFreq <- ifelse(comm.sum[,] >= 0,1,0) ## turning into a presense (1) absence (0) dataframe

PlotRichness <- as.data.frame(apply(SpeciesFreq, 1, sum, na.rm = TRUE)) ## recalculating a richness per plot
colnames(PlotRichness) <- "richness"

SpeciesFreq <- as.data.frame(apply(SpeciesFreq, 2, sum, na.rm = TRUE)) ## getting the number of times a species occurred 
colnames(SpeciesFreq) <- "SpeciesFrequency" ## renaming to be be easy to identify later
SpeciesFreq$Species <- rownames(SpeciesFreq)

rm(tmp);rm(i);rm(vec);rm(plot.avg);rm(trans.avg)

### Should sites be split in analysis? ####
## PERMANOVA 
comm.sum[comm.sum == 0] <- 0.0001
comm.sum[is.na(comm.sum)] <- 0

comm.dat <- as.matrix(comm.sum) ## making a sep df to test comm. data
soils <- soils[order(soils$Field.ID),]
soils.dat <- as.matrix(soils[,c(5:16)])
rownames(soils.dat) <- soils$Field.ID
env <- env[order(env$plot),]
sites <- factor(env$site) ## the sites are not perfectly ordered so I have to use these numbers
sev <- factor(env$sev, levels = c("unburn","low","mod","high"))
trt <- factor(env$trt)

PERMANOVA.comm <- adonis2(comm.dat ~ sites+sev+trt, 
                          permutations = 999,
                          method = "bray")
PERMANOVA.soil <- adonis2(soils.dat ~ sites+sev+trt, 
                          permutations = 999,
                          method = "bray")

print(PERMANOVA.comm) ## differences between sites
print(PERMANOVA.soil) ## differences between sites

# Check homogeneity of dispersions (important assumption)
# This tests if variance within groups is similar
dist_matrix.comm <- vegdist(comm.dat, method = "bray")
dist_matrix.soil <- vegdist(soils.dat, method = "bray")

dispersion.comm <- betadisper(dist_matrix.comm, sites)
dispersion.soil <- betadisper(dist_matrix.soil, sites)

permutest(dispersion.comm, permutations = 999)
permutest(dispersion.soil, permutations = 999)

# Visualize dispersions
par(mfrow = c(1,2))
plot(dispersion.comm) ## sites are very different
plot(dispersion.soil) ## sites are very different

rm(PERMANOVA.comm);rm(dist_matrix.comm);rm(dispersion.comm)
rm(PERMANOVA.soil);rm(dist_matrix.soil);rm(dispersion.soil)
rm(comm.dat);rm(soils.dat);rm(trt);rm(sev);rm(sites)

#### Question 1 ####
## SOILS Data question
## How does fire retardant application influence soil nitrogen and 
## phosphorous compared to unamended burned and unburned areas?

## preliminary visualizations
## looking at variation between sites, treatments, and burn severity

par(mfrow = c(1,1))
## NH4
plot(soils$NH4, col = as.factor(soils$Location), pch = 16) 
plot(soils$NH4, col = as.factor(soils$Treatment), pch = 16)
plot(soils$NH4, col = as.factor(soils$Burn.Severity), pch = 16)
## seems more related to burn severity

## NO3
plot(soils$NO3, col = as.factor(soils$Location), pch = 16)
plot(soils$NO3, col = as.factor(soils$Treatment), pch = 16)
plot(soils$NO3, col = as.factor(soils$Burn.Severity), pch = 16)
## NO3 is higher in trt regardless of burn severity

## PO4
plot(soils$PO4, col = as.factor(soils$Location), pch = 16)
plot(soils$PO4, col = as.factor(soils$Treatment), pch = 16)
plot(soils$PO4, col = as.factor(soils$Burn.Severity), pch = 16)
## looks like heavily related to trt

## splitting into two different groups
quarry.soils <- soils[soils$Location == "Quarry",]
quarry.soils$Burn.Severity <- factor(quarry.soils$Burn.Severity, levels = c("Unburned", "Low", "Moderate", "High"))
quarry.soils$Treatment <- as.factor(quarry.soils$Treatment)

sc.soils <- soils[soils$Location == "Stone Canyon",]
sc.soils$Burn.Severity[sc.soils$Burn.Severity == "Fire  Line"] <- "Burned"
sc.soils$Burn.Severity <- factor(sc.soils$Burn.Severity, levels = c("Unburned", "Burned"))
sc.soils$Treatment <- as.factor(sc.soils$Treatment)

## Linear Models
## NH4
NH4.mod.q <- lm(NH4 ~ Burn.Severity*Treatment, data = quarry.soils)
summary(NH4.mod.q)
plot(residuals(NH4.mod.q))
qqnorm(residuals(NH4.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NH4.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

NH4.mod.sc <- lm(NH4 ~ Burn.Severity*Treatment, data = sc.soils)
summary(NH4.mod.sc)
plot(residuals(NH4.mod.sc))
qqnorm(residuals(NH4.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NH4.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough
rm(NH4.mod.q);rm(NH4.mod.sc)

## NO3
NO3.mod.q <- lm(NO3 ~ Burn.Severity*Treatment, data = quarry.soils)
summary(NO3.mod.q)
plot(residuals(NO3.mod.q))
qqnorm(residuals(NO3.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NO3.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

NO3.mod.sc <- lm(NO3 ~ Burn.Severity*Treatment, data = sc.soils)
summary(NO3.mod.sc)
plot(residuals(NO3.mod.sc))
qqnorm(residuals(NO3.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NO3.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough
rm(NO3.mod.q);rm(NO3.mod.sc)

## PO4
P.mod.q <- lm(PO4 ~ Burn.Severity*Treatment, data = quarry.soils)
summary(P.mod.q)
plot(residuals(P.mod.q))
qqnorm(residuals(P.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(P.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

P.mod.sc <- lm(PO4 ~ Burn.Severity*Treatment, data = sc.soils)
summary(P.mod.sc)
plot(residuals(P.mod.sc))
qqnorm(residuals(P.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(P.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough
rm(P.mod.q);rm(P.mod.sc)

## supplemental analysis of additional nutrients
colnames(soils[c(5:16)])

## pH
pH.mod.q <- lm(pH ~ Burn.Severity*Treatment, data = quarry.soils)
summary(pH.mod.q)
pH.mod.sc <- lm(pH ~ Burn.Severity*Treatment, data = sc.soils)
summary(pH.mod.sc)
rm(pH.mod.q);rm(pH.mod.sc)

## Na
Na.mod.q <- lm(Na ~ Burn.Severity*Treatment, data = quarry.soils)
summary(Na.mod.q)
Na.mod.sc <- lm(Na ~ Burn.Severity*Treatment, data = sc.soils)
summary(Na.mod.sc)
rm(Na.mod.q);rm(Na.mod.sc)

## K
K.mod.q <- lm(K ~ Burn.Severity*Treatment, data = quarry.soils)
summary(K.mod.q)
K.mod.sc <- lm(K ~ Burn.Severity*Treatment, data = sc.soils)
summary(K.mod.sc)
rm(K.mod.q);rm(K.mod.sc)

## SO4
SO4.mod.q <- lm(SO4 ~ Burn.Severity*Treatment, data = quarry.soils)
summary(SO4.mod.q)
SO4.mod.sc <- lm(SO4 ~ Burn.Severity*Treatment, data = sc.soils)
summary(SO4.mod.sc)
rm(SO4.mod.q);rm(SO4.mod.sc)

## Looking into how much the nutrients vary across treatments
aggregate(quarry.soils$NH4 ~ quarry.soils$Treatment, FUN= mean)
aggregate(sc.soils$NH4 ~ sc.soils$Treatment, FUN= mean)

aggregate(quarry.soils$NO3 ~ quarry.soils$Treatment, FUN= mean) ## order of magnitude greater
aggregate(sc.soils$NO3 ~ sc.soils$Treatment, FUN= mean)

aggregate(quarry.soils$PO4 ~ quarry.soils$Treatment, FUN= mean) ## order of magnitude greater
aggregate(sc.soils$PO4 ~ sc.soils$Treatment, FUN= mean) ## two orders of magnitude greater
## consistently higher soil nutrient concentrations in the fire retardant treatments

#### Figure 1 - Soils ####
soils$Burn.Severity[soils$Burn.Severity == "Fire  Line"] <- "Burned"
soils$Burn.Severity[soils$Burn.Severity == "Burned"] <- "Low"

se <- function(x){sd(x)/sqrt(length(x))} ## creating a function for standard error
soils$plotting <- ifelse(soils$Burn.Severity == "Unburned" & soils$Treatment == "Control", "UCON",
                                ifelse(soils$Burn.Severity == "Unburned" & soils$Treatment == "Fire Retardant", "UFR",
                                       ifelse(soils$Burn.Severity == "Low" & soils$Treatment == "Control", "LCON",
                                              ifelse(soils$Burn.Severity == "Low" & soils$Treatment == "Fire Retardant", "LFR", 
                                                     ifelse(soils$Burn.Severity == "Moderate" & soils$Treatment == "Control", "MCON",
                                                            ifelse(soils$Burn.Severity == "Moderate" & soils$Treatment == "Fire Retardant", "MFR",
                                                                   ifelse(soils$Burn.Severity == "High" & soils$Treatment == "Control", "HCON",
                                                                          ifelse(soils$Burn.Severity == "High" & soils$Treatment == "Fire Retardant", "HFR", NA))))))))
## setting figure parameters 
Fig2order.sc <- c("UCON","LCON")
vec1.sc <- c(0.9,1.9)
Fig2paired.sc <- c("UFR","LFR")
vec2.sc <- c(1.1,2.1)

Fig2order.q <- c("UCON","LCON","MCON","HCON")
vec1.q <- c(0.9,1.9,2.9,3.9)
Fig2paired.q <- c("UFR","LFR","MFR","HFR")
vec2.q <- c(1.1,2.1,3.1,4.1)

SC <- soils[soils$Location == "Stone Canyon",]
Quarry <- soils[soils$Location == "Quarry",]

par(mfrow = c(3,2))

## change the nutrients of interest based on the column name in soils df
nut <- 'NH4'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'NO3'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'PO4'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

#### Supplemental Soil Nutrient Figures ####
## pH
## change the nutrients of interest based on the column name in soils df
par(mfrow = c(1,2))
nut <- 'pH'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 7, lty = 2)
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 7, lty = 2)
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## Na, K, SO4
par(mfrow = c(3,2))

## change the nutrients of interest based on the column name in soils df
nut <- 'Na'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'K'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'SO4'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(soils[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

rm(Quarry);rm(SC);rm(Fig2order.q);rm(Fig2paired.q);rm(Fig2order.sc);rm(Fig2paired.sc)
rm(i);rm(vec1.sc);rm(vec2.sc);rm(vec1.q);rm(vec2.q)
rm(nut)
par(mfrow = c(1,1))


#### Question 2  ####
## Community Composition Question
## How do vegetative communities in fire-retardant drop zones differ in species richness, 
## community composition, and cover compared to burned, unburned, and adjacent fire 
## perimeter sites without drops?

## part 1 - general trend cover/richness (boxplots)  - supplement split them
## part 2 - general trend cover/richness (native vs. invasive spp.) - supplement split them
## part 3 - glmm broken out by site (look at site differences)
## part 4 - ordination of whole comm (site being a split)

## part 1.1 - richness between control and trt
PlotRichness$sev <- env$sev[match(rownames(PlotRichness),env$plot)]
PlotRichness$trt <- env$trt[match(rownames(PlotRichness),env$plot)]
PlotRichness$site <- env$site[match(rownames(PlotRichness),env$plot)]

str(PlotRichness)
t.test(PlotRichness$richness~ PlotRichness$site, na.rm = TRUE)
## difference between the two sites (quarry has fewer species on average)

## part 1.2  - cover between control and trt
PlotCov$sev <- env$sev[match(rownames(PlotCov),env$plot)]
PlotCov$trt <- env$trt[match(rownames(PlotCov),env$plot)]
PlotCov$site <- env$site[match(rownames(PlotCov),env$plot)]

str(PlotCov)
t.test(PlotCov$TotalPlotCover~ PlotCov$site, na.rm = TRUE)
## no difference in plot cover between sites

## part 2 - native and invasive cover
## separate out native vs. invasive species for comparison across burn severities (basically do part 1 but 4 times)
## plot level proportion richness Invasive /proportion cover invasive
## functional group and life history can be added if interested

status <- comm.sum
colnames(status)
match(colnames(status), sp.info$code)

native <- status[,colnames(status) %in% sp.info$code[sp.info$status == "N"]]
invasive <- status[,colnames(status) %in% sp.info$code[sp.info$status == "I"]]
annual <- status[,colnames(status) %in% sp.info$code[sp.info$duration == "annual"]]
perennial <- status[,colnames(status) %in% sp.info$code[sp.info$duration == "perennial"]]

env$X <- NULL

native.cov <- data.frame(native.cov = apply(native, 1 , sum),
                         plot = rownames(native))
invasive.cov <- data.frame(invasive.cov = apply(invasive,1,sum),
                           plot = rownames(invasive))
annual.cov <- data.frame(annual.cov = apply(annual, 1 , sum),
                         plot = rownames(annual))
perennial.cov <- data.frame(perennial.cov = apply(perennial, 1 , sum),
                         plot = rownames(perennial))
env$native.cov <- native.cov$native.cov[match(native.cov$plot, env$plot)]
env$invasive.cov <- invasive.cov$invasive.cov[match(invasive.cov$plot, env$plot)]
env$annual.cov <- annual.cov$annual.cov[match(annual.cov$plot, env$plot)]
env$perennial.cov <- perennial.cov$perennial.cov[match(perennial.cov$plot, env$plot)]

env$prop.inv.cov <- env$invasive.cov/(env$native.cov+env$invasive.cov)
env$prop.ann.cov <- env$annual.cov/(env$perennial.cov+env$annual.cov)
hist(env$prop.inv.cov)
hist(env$prop.ann.cov)

grad <- native.cov
grad$invasive.cov <- invasive.cov$invasive.cov[match(grad$plot, invasive.cov$plot)]
grad$annual.cov <- annual.cov$annual.cov[match(grad$plot, annual.cov$plot)]
grad$plot.cov <- PlotCov$TotalPlotCover[match(grad$plot, rownames(PlotCov))]
grad$plot.rich <- PlotRichness$richness[match(grad$plot, rownames(PlotRichness))]
grad$NO3 <- soils$NO3[match(grad$plot, soils$Field.ID)]
grad$NH4 <- soils$NH4[match(grad$plot, soils$Field.ID)]
grad$PO4 <- soils$PO4[match(grad$plot, soils$Field.ID)]
grad$site <- soils$Location[match(grad$plot, soils$Field.ID)]
grad$prop.inv <- grad$invasive.cov/grad$plot.cov
grad$prop.ann <- grad$annual.cov/grad$plot.cov
grad$trt <- PlotCov$trt[match(grad$plot, rownames(PlotCov))]
grad$sev <- PlotCov$sev[match(grad$plot, rownames(PlotCov))]

rm(native.cov);rm(invasive.cov)

native <- ifelse(native > 0, 1,0)
invasive <- ifelse(invasive > 0,1,0)
annual <- ifelse(annual > 0, 1,0)
perennial <- ifelse(perennial > 0,1,0)

env$native.sp <- apply(native, 1 , sum)
env$invasive.sp <- apply(invasive, 1, sum)
env$annual.sp <- apply(annual, 1 , sum)
env$perennial.sp <- apply(perennial, 1, sum)

env$prop.inv.sp <- env$invasive.sp/(env$native.sp+env$invasive.sp)
grad$prop.inv.sp <- env$prop.inv.sp[match(grad$plot, env$plot)]
env$prop.ann.sp <- env$annual.sp/(env$perennial.sp+env$annual.sp)
grad$prop.ann.sp <- env$prop.ann.sp[match(grad$plot, env$plot)]
hist(env$prop.inv.sp)
hist(env$prop.ann.sp)

hist(env$prop.inv.sp[env$trt == "con"])
hist(env$prop.inv.sp[env$trt == "fr"])

hist(env$prop.ann.sp[env$trt == "con"])
hist(env$prop.ann.sp[env$trt == "fr"])

#### Linear Models to address Q2 ####
grad$trt <- factor(grad$trt, levels = c("con", "fr"))
grad$sev <- factor(grad$sev, levels = c("unburn", "low", "mod", "high"))
str(grad)
pseudo.R.squared <- function(model){
  1 - (model$deviance/model$null.deviance)
}
q.grad <- grad[grad$site == "Quarry",]
sc.grad <- grad[grad$site == "Stone Canyon",]

## Quarry models
## response ~ soils + trt + severity
summary(lm(plot.cov ~ NO3+NH4+PO4+trt+sev,data = q.grad)) ## keeping in the trt and sev to control for it
summary(lm(prop.inv ~ NO3+NH4+PO4+trt+sev,data = q.grad))
summary(lm(prop.ann ~ NO3+NH4+PO4+trt+sev,data = q.grad))
summary(glm(plot.rich ~ NO3+NH4+PO4+trt+sev,data = q.grad, family = poisson(link = "log")))
pseudo.R.squared(glm(plot.rich ~ NO3+NH4+PO4+trt+sev,data = q.grad, family = poisson(link = "log")))
summary(lm(prop.inv.sp ~ NO3+NH4+PO4+trt+sev,data = q.grad))
summary(lm(prop.ann.sp ~ NO3+NH4+PO4+trt+sev,data = q.grad))


## Stone canyon
## response ~ soils + trt + severity
summary(lm(plot.cov ~ NO3+NH4+PO4+trt+sev,data = sc.grad)) ## keeping in the trt and sev to control for it
summary(lm(prop.inv ~ NO3+NH4+PO4+trt+sev,data = sc.grad))
summary(lm(prop.ann ~ NO3+NH4+PO4+trt+sev,data = sc.grad))
summary(glm(plot.rich ~ NO3+NH4+PO4+trt+sev,data = sc.grad, family = poisson(link = "log")))
pseudo.R.squared(glm(plot.rich ~ NO3+NH4+PO4+trt+sev,data = sc.grad, family = poisson(link = "log")))
summary(lm(prop.inv.sp ~ NO3+NH4+PO4+trt+sev,data = sc.grad))
summary(lm(prop.ann.sp ~ NO3+NH4+PO4+trt+sev,data = sc.grad))


#### Figure 2 - Richness ####
se <- function(x){sd(x)/sqrt(length(x))} ## creating a function for standard error
sc.grad$plotting <- ifelse(sc.grad$sev == "unburn" & sc.grad$tr == "con", "UCON",
                           ifelse(sc.grad$sev == "unburn" & sc.grad$tr == "fr", "UFR",
                                  ifelse(sc.grad$sev == "low" & sc.grad$tr == "con", "LCON",
                                         ifelse(sc.grad$sev == "low" & sc.grad$tr == "fr", "LFR", 
                                                ifelse(sc.grad$sev == "mod" & sc.grad$tr == "con", "MCON",
                                                       ifelse(sc.grad$sev == "mod" & sc.grad$tr == "fr", "MFR",
                                                              ifelse(sc.grad$sev == "high" & sc.grad$tr == "con", "HCON",
                                                                     ifelse(sc.grad$sev == "high" & sc.grad$tr == "fr", "HFR", NA))))))))
q.grad$plotting <- ifelse(q.grad$sev == "unburn" & q.grad$tr == "con", "UCON",
                           ifelse(q.grad$sev == "unburn" & q.grad$tr == "fr", "UFR",
                                  ifelse(q.grad$sev == "low" & q.grad$tr == "con", "LCON",
                                         ifelse(q.grad$sev == "low" & q.grad$tr == "fr", "LFR", 
                                                ifelse(q.grad$sev == "mod" & q.grad$tr == "con", "MCON",
                                                       ifelse(q.grad$sev == "mod" & q.grad$tr == "fr", "MFR",
                                                              ifelse(q.grad$sev == "high" & q.grad$tr == "con", "HCON",
                                                                     ifelse(q.grad$sev == "high" & q.grad$tr == "fr", "HFR", NA))))))))

## setting figure parameters 
Fig2order.sc <- c("UCON","LCON")
vec1.sc <- c(0.9,1.9)
Fig2paired.sc <- c("UFR","LFR")
vec2.sc <- c(1.1,2.1)

Fig2order.q <- c("UCON","LCON","MCON","HCON")
vec1.q <- c(0.9,1.9,2.9,3.9)
Fig2paired.q <- c("UFR","LFR","MFR","HFR")
vec2.q <- c(1.1,2.1,3.1,4.1)

SC <- sc.grad
Quarry <- q.grad

par(mfrow = c(3,2))

## change the nut variable (borrowed from soils code) based on the column name in grad df
nut <- 'plot.rich'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(grad[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(grad[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'prop.inv.sp'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
abline(h = 0.5, lty = 2)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'prop.ann.sp'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

#### Figure 3 - Cover ####
par(mfrow = c(3,2))

## change the nut variable (borrowed from soils code) based on the column name in grad df
nut <- 'plot.cov'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,max(grad[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,max(grad[,nut])),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'prop.inv'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
abline(h = 0.5, lty = 2)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

## change the nutrients of interest based on the column name in soils df
nut <- 'prop.ann'
plot(x = c(0.5:2.5),
     y = c(0.5:2.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:2), line = 1, tick = F, labels = c("Unburned", "Low"), cex.axis = 1.5)
for(i in 1:length(Fig2order.sc)){
  points(x = rep(vec1.sc[i], length(SC[,nut][SC$plotting == Fig2order.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2order.sc[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.sc[i], length(SC[,nut][SC$plotting == Fig2paired.sc[i]])),
         y = SC[,nut][SC$plotting == Fig2paired.sc[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2order.sc[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])-se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x0 = vec1.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2order.sc[i]])+se(SC[,nut][SC$plotting == Fig2order.sc[i]])), x1 = vec1.sc[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.sc[i]+0.1,
         y = mean(SC[,nut][SC$plotting == Fig2paired.sc[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])-se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x0 = vec2.sc[i]+0.1, 
           y1 = (mean(SC[,nut][SC$plotting == Fig2paired.sc[i]])+se(SC[,nut][SC$plotting == Fig2paired.sc[i]])), x1 = vec2.sc[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

plot(x = c(0.5:4.5),
     y = c(0.5:4.5),
     ylim = c(0,1),
     las = 1,
     cex.axis = 1.5,
     ylab = "", ## density
     type = "n",
     xaxt = "n",
     xlab = "") ## disturbance history
abline(h = 0.5, lty = 2)
axis(1, at = c(1:4), line = 1, tick = F, labels = c("Unburned", "Low", "Moderate", "High"), cex.axis = 1.5)
for(i in 1:length(Fig2order.q)){
  points(x = rep(vec1.q[i], length(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2order.q[i]],
         col = rgb(0,0,0, alpha = 0.5),
         pch = 19)
  points(x = rep(vec2.q[i], length(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])),
         y = Quarry[,nut][Quarry$plotting == Fig2paired.q[i]],
         col = rgb(1,0,0, alpha = 0.5),
         pch = 19)
  
  points(x = vec1.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]]),
         col = rgb(0,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x0 = vec1.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2order.q[i]])), x1 = vec1.q[i]+0.1, 
           col = rgb(0,0,0),lwd = 1.5)
  
  points(x = vec2.q[i]+0.1,
         y = mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]]),
         col = rgb(1,0,0, alpha = 0.5),
         pch = 17)
  segments(y0 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])-se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x0 = vec2.q[i]+0.1, 
           y1 = (mean(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])+se(Quarry[,nut][Quarry$plotting == Fig2paired.q[i]])), x1 = vec2.q[i]+0.1, 
           col = rgb(1,0,0),lwd = 1.5)
  
}

rm(Quarry);rm(SC);rm(Fig2order.q);rm(Fig2paired.q);rm(Fig2order.sc);rm(Fig2paired.sc)
rm(i);rm(vec1.sc);rm(vec2.sc);rm(vec1.q);rm(vec2.q)
rm(nut)


# #### Supplemental - cov ~ PO4 ####
# par(mfrow = c(1,1))
# 
# mod <- lm(plot.cov ~ NO3+NH4+PO4+trt+sev, data = q.grad)
# plot(plot.cov ~ PO4, data = q.grad,
#      las = 1,
#      log = "y",
#      cex.axis = 1.5,
#      pch = 16,
#      cex = 0.75,
#      xlab = " ",
#      ylab = " ")
# newdata <- data.frame(NO3=mean(q.grad$NO3, na.rm = T),
#                       NH4=mean(q.grad$NH4, na.rm = T),
#                       PO4=seq(min(q.grad$PO4, na.rm = T), max(q.grad$PO4, na.rm = T), length.out = 30),
#                       trt=q.grad$trt[1],
#                       sev=q.grad$sev[1])
# preds <- predict(mod, newdata, type = "response", se.fit = TRUE)
# preds$upperci <- preds$fit + (1.96*preds$se.fit)
# preds$lowerci <- preds$fit - (1.96*preds$se.fit)
# polygon(x = c(newdata$PO4, rev(newdata$PO4)),
#         y = c(preds$lowerci, rev(preds$upperci)),
#         col = adjustcolor("black", alpha.f = 0.25),
#         border = NA)
# lines(newdata$PO4, preds$fit, lty = 1)
# 

#### Question 2 - Supplemental Ordinations ####
str(env)
env$sev <- factor(env$sev, levels = c("unburn", "low", "mod", "high"))
env$trt <- as.factor(env$trt)
levels(env$sev)
env$NO3 <- soils$NO3[match(env$plot, soils$Field.ID)]
env$NH4 <- soils$NH4[match(env$plot, soils$Field.ID)]
env$PO4 <- soils$PO4[match(env$plot, soils$Field.ID)]
env$tot.cov <- PlotCov$TotalPlotCover[match(env$plot, rownames(PlotCov))]
env$rich <- PlotRichness$richness[match(env$plot, rownames(PlotRichness))]


## Quarry Fire  
Quarry <- as.matrix(comm.sum[match(env$plot[env$site == "quarry"], rownames(comm.sum)),])
str(Quarry)
env.q <- env[match(rownames(Quarry),env$plot),]

## Stone Canyon 
SC <- as.matrix(comm.sum[match(env$plot[env$site == "sc"], rownames(comm.sum)),])
str(SC)
env.sc <- env[match(rownames(SC),env$plot),]

pal <- c("#f6d746", "#e55c30", "#84206b", "#140b34")
pal2 <- c("black", "red")


## Quarry NMDS Ordination
ord.nmds.stress.q <- rep(NA,10)
for(i in 1:10){
  ord.nmds.stress.q[i] <- metaMDS(Quarry, k = i, try = 1000, distance = "bray")$stress 
}
par(mfrow = c(2,2))
plot(ord.nmds.stress.q,
     ylim = c(0,1),
     ylab = "Stress",
     xlab = "Number of Dimensions",
     main = "Quarry Fire",
     las = 1,
     pch = 16)
abline(h = 0.2)

set.seed(1)
NMDSord.q <- metaMDS(Quarry, k = 3, try = 1000, distance = "bray") ## convergence
stressplot(NMDSord.q,
           main = "Quarry Fire")
## Stone Canyon NMDS Ordination
ord.nmds.stress.sc <- rep(NA,10)
for(i in 1:10){
  ord.nmds.stress.sc[i] <- metaMDS(SC, k = i, try = 1000, distance = "bray")$stress 
}
plot(ord.nmds.stress.sc,
     ylim = c(0,1),
     ylab = "Stress",
     xlab = "Number of Dimensions",
     main = "Stone Canyon",
     las = 1,
     pch = 16)
abline(h = 0.2)

set.seed(1)
NMDSord.sc <- metaMDS(SC, k = 4, try = 1000, distance = "bray") ## convergence
stressplot(NMDSord.sc,
           main = "Stone Canyon")

par(mfrow = c(2,2))
## a - severity
plot(NMDSord.q, type = "n",display = "sites", las = 1) 
points(NMDSord.q, display = "sites", pch = 19, cex = .75, col = pal[as.factor(env.q$sev)])
vectors <- envfit(NMDSord.q, env.q[,c(11,14:19)],na.rm = TRUE)
ordiellipse(NMDSord.q, display = "sites", env.q$sev, draw = "lines",
            col = pal, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Unburned", "Low", "Mod", "High"), col = pal, pch = 15, cex = 1, ncol = 2, bty = "n")

## b - fire retardant
plot(NMDSord.q, type = "n",display = "sites", las = 1) 
points(NMDSord.q, display = "sites", pch = 19, cex = .75, col = pal2[as.factor(env.q$trt)])
vectors <- envfit(NMDSord.q, env.q[,c(11,14:19)],na.rm = TRUE)
ordiellipse(NMDSord.q, display = "sites", env.q$trt, draw = "lines",
            col = pal2, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Control", "Fire Retardant"), col = pal2, pch = 15, cex = 1, ncol = 2, bty = "n")

## c - severity
plot(NMDSord.sc, type = "n",display = "sites", las = 1) 
points(NMDSord.sc, display = "sites", pch = 19, cex = .75, col = pal[as.factor(env.sc$sev)])
vectors <- envfit(NMDSord.sc, env.sc[,c(11,14:19)],na.rm = TRUE)
ordiellipse(NMDSord.sc, display = "sites", env.sc$sev, draw = "lines",
            col = pal, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Unburned", "Low", "Mod", "High"), col = pal, pch = 15, cex = 1, ncol = 2, bty = "n")

## d - fire retardant
plot(NMDSord.sc, type = "n",display = "sites", las = 1) 
points(NMDSord.sc, display = "sites", pch = 19, cex = .75, col = pal2[as.factor(env.sc$trt)])
vectors <- envfit(NMDSord.sc, env.sc[,c(11,14:19)],na.rm = TRUE)
ordiellipse(NMDSord.sc, display = "sites", env.sc$trt, draw = "lines",
            col = pal2, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Control", "Fire Retardant"), col = pal2, pch = 15, cex = 1, ncol = 2, bty = "n")

rm(i);rm(ord.nmds.stress.q);rm(pal);rm(pal2);rm(se);rm(vectors);rm(Quarry)
rm(SC);rm(env.q);rm(env.sc);rm(NMDSord.q)
rm(NMDSord.sc);rm(ord.nmds.stress.sc);rm(mod);rm(preds);rm(newdata)