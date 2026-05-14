#### SET UP #####
##### Import Libraries #####
library(dplyr) ## useful for transforming data
library(stringr) ## needed to find patterns in character strings
library(vegan) ## multidimentional vegetation analysis package
library(piecewiseSEM) # for structural equation modeling
library(semEffect) ## get direct and indirect effects from piecewise SEM
library(MASS) # for negative binomial glm
library(modelbased) 
library(patchwork) # for combining ggplots
library(here) # for file management
library(ggplot2) # for plotting
library(magrittr) # for the pipe operator
library(tidyverse) # for data manipulation and visualization
library(car) # for vi function to check for multicollinearity
library(finalfit) #round tidy

##### Set up project structure #####
setwd(here()) ## set working directory 
dir.create(here("output")) ## create an output folder for figures and tables)

##### Set ggplot2 custom plotting theme
theme_new <- function(base_size = 9,base_family = "helvetica"){
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line.x = element_line(color="black", linewidth = 0.25),
      axis.line.y = element_line(color="black", linewidth = 0.25),
      axis.title = element_text(size = 9),
      axis.text = element_text(colour="black", size=8),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.grid = element_blank(),   
      plot.background = element_rect(fill = NA, colour = NA),
      panel.border = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 9)
    )
}
theme_set(theme_new())

#### DATA IMPORT, CLEANING, AND SUMMARIZATION #####
##### Reading in Data #####
cover <- read.csv("./data/FieldData(comm).csv") ## field data on plant cover
rich <- read.csv("./data/FieldData(rich).csv") ## field data on species richness
env <- read.csv("./data/ENV_Data.csv") ## environmental data from script 1-Clean Environmental Data.R
sp.info <- read.csv("./data/IntroducedStatusKey.csv") ## species information on introduced/native, functional group, etc. 
soils <- read.csv("./data/CarterFR_analysis.csv")

##### Cleaning Data #####
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

##### Calculate Average Species Cover per Plot #####
table(cover$plot) ## looking at plots and double checking
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

##### Creating Community Data (Cover Estimates + Rare Species only found as richness) #####
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


#### INITIAL EXPLORATION OF THE DATA #####
##### Should sites be split in analysis? ######
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
jpeg(here("output/DispersionPlots.jpeg"), width = 700, height = 400)
par(mfrow = c(1,2),cex=1.2)
plot(dispersion.comm, col=c("#018571","#a6611a"), main="plant community") ## sites are very different
plot(dispersion.soil, col=c("#018571","#a6611a"), main="soils") ## sites are very different
dev.off()

rm(PERMANOVA.comm);rm(dist_matrix.comm);rm(dispersion.comm)
rm(PERMANOVA.soil);rm(dist_matrix.soil);rm(dispersion.soil)
rm(comm.dat);rm(soils.dat);rm(trt);rm(sev);rm(sites)


# How does fire retardant application influence soil nitrogen and 
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


#### ANALYSES ####
##### Combine data ####
dat.cover <- comm.sum.long %>% 
  left_join(., sp.info, by = c("code" = "code")) %>% 
  group_by(plot) %>% 
  summarize(plot.cov = sum(plotavg, na.rm=T),
            invasive.cov = sum(plotavg[status == "I"], na.rm=T),
            annual.cov = sum(plotavg[duration== "annual"], na.rm=T)) %>% 
  mutate(prop.inv.cov = invasive.cov/plot.cov,
         prop.ann.cov = annual.cov/plot.cov)

dat.rich <- comm.sum.long %>% 
  left_join(., sp.info, by = c("code" = "code")) %>% 
  group_by(plot) %>% 
  summarize(plot.rich = n(),
            plot.notCarex = sum(code != "CAREX SP.", na.rm=T),
            invasive.rich = sum(status=="I", na.rm=T), 
            annual.rich = sum(duration=="annual", na.rm=T)) %>%
  mutate(prop.inv.sp = invasive.rich/plot.notCarex,
         prop.ann.sp = annual.rich/plot.notCarex) %>% 
  dplyr::select(-plot.notCarex)

dat <- env[,c("plot", "site", "sev", "trt")] %>% 
  full_join(., soils[,c(1,5:16)], by = c("plot" = "Field.ID")) %>% 
  full_join(., dat.cover, by = "plot") %>% 
  full_join(., dat.rich, by = "plot") %>% 
  filter(!is.na(NH4)) %>% 
  # change treatment to binary and burn severity to numeric
  mutate(Treatment=ifelse(trt=="con", 0, 1), 
         Burn.Severity=as.numeric(factor(sev, levels=c("unburn", "low", "mod", "high"), ordered=T)))

##### QUARRY #####
q.dat <- dat %>% 
  filter(site=="quarry") 
rownames(q.dat) <- 1:nrow(q.dat) # this needs to be run in order for the standardized coefficients to be calculated

###### Soils ######
## NH4
NH4.mod.q <- lm(NH4 ~ Burn.Severity*Treatment, data = q.dat)
summary(NH4.mod.q)
plot(residuals(NH4.mod.q))
qqnorm(residuals(NH4.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NH4.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## NO3
NO3.mod.q <- lm(NO3 ~ Burn.Severity*Treatment, data =  q.dat)
summary(NO3.mod.q)
plot(residuals(NO3.mod.q))
qqnorm(residuals(NO3.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NO3.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## PO4
P.mod.q <- lm(PO4 ~ Burn.Severity*Treatment, data =  q.dat)
summary(P.mod.q)
plot(residuals(P.mod.q))
qqnorm(residuals(P.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(P.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

###### COVER ######
## Total
cov.mod.q  <- lm(plot.cov ~  NO3 + NH4 + PO4 + Burn.Severity + Treatment, data = q.dat)
summary(cov.mod.q)
vif(cov.mod.q)
plot(residuals(cov.mod.q))
qqnorm(residuals(cov.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(cov.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion invasive 
icov.mod.q  <- lm(prop.inv.cov  ~  NO3 + NH4 + PO4 + Burn.Severity + Treatment, data = q.dat)
summary(icov.mod.q )
plot(residuals(icov.mod.q ))
qqnorm(residuals(icov.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(icov.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion annual
acov.mod.q  <- lm(prop.ann.cov ~  NO3 + NH4 + PO4 + Burn.Severity + Treatment, data = q.dat)
summary(acov.mod.q )
plot(residuals(acov.mod.q ))
qqnorm(residuals(acov.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(acov.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

###### RICHNESS #####
## Total richness
rich.mod.q <- glm.nb(plot.rich ~  NO3 + NH4 + PO4 +  Burn.Severity + Treatment, data = q.dat) # mean(dat$plot.rich) = 20; var(dat$plot.rich) = 54 
summary(rich.mod.q)
plot(residuals(rich.mod.q ))
qqnorm(residuals(rich.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(rich.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion invasive richnes
irich.mod.q <- lm(prop.inv.sp~  NO3 + NH4 + PO4 + Burn.Severity + Treatment, data = q.dat) 
summary(irich.mod.q)
plot(residuals(irich.mod.q ))
qqnorm(residuals(irich.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(irich.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion annual richness
arich.mod.q <- lm(prop.ann.sp~  NO3 + NH4 + PO4 + Burn.Severity + Treatment, data = q.dat) 
summary(arich.mod.q)
plot(residuals(arich.mod.q ))
qqnorm(residuals(arich.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(arich.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

###### PIECEWISE SEM BY PLANT RESPONSE #####
## total cover
modList.plot.cov <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, cov.mod.q, data =q.dat)
summary(modList.plot.cov)
coef.summary <- summary(modList.plot.cov)$coefficients
di.effect <- get_effect(modList.plot.cov, target="plot.cov")$effect_long %>% mutate(Response="total cover", site="Quarry")
r2.summary <- summary(modList.plot.cov)$R2 %>% mutate(site="Quarry") 
Cstat.summary <- summary(modList.plot.cov)$Cstat  %>% mutate(Response="total cover", site="Quarry")

## invasive cover
modList.plot.icov <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, icov.mod.q, data =q.dat)
summary(modList.plot.icov)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.icov)$coefficients) 
di.effect <- bind_rows(di.effect, get_effect(modList.plot.icov, target="prop.inv.cov")$effect_long %>% mutate(Response="prop.invasive cover", site="Quarry"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.icov)$R2 %>% mutate(site="Quarry")) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.icov)$Cstat  %>% mutate(Response="prop. invasive cover", site="Quarry"))

## annual cover
modList.plot.acov <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, acov.mod.q, data =q.dat)
summary(modList.plot.acov)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.acov)$coefficients) 
di.effect <- bind_rows(di.effect, get_effect(modList.plot.acov, target="prop.ann.cov")$effect_long %>% mutate(Response="prop. annual cover", site="Quarry"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.acov)$R2 %>% mutate(site="Quarry") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.acov)$Cstat  %>% mutate(Response="prop. annual cover", site="Quarry"))

## total richness
modList.plot.rich <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, rich.mod.q, data =q.dat)
summary(modList.plot.rich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.rich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.rich, target="plot.rich")$effect_long %>% mutate(Response="species richness", site="Quarry"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.rich)$R2 %>% mutate(site="Quarry") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.rich)$Cstat  %>% mutate(Response="species richness", site="Quarry"))

## annual rich
modList.plot.arich <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, arich.mod.q, data =q.dat)
summary(modList.plot.arich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.arich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.arich, target="prop.ann.sp")$effect_long %>% mutate(Response="prop. annual richness", site="Quarry"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.arich)$R2 %>% mutate(site="Quarry") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.arich)$Cstat  %>% mutate(Response="prop. annual richness", site="Quarry"))

## invasive rich
modList.plot.irich <- psem(NH4.mod.q, NO3.mod.q, P.mod.q, irich.mod.q, data =q.dat)
summary(modList.plot.irich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.irich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.irich, target="prop.inv.sp")$effect_long %>% mutate(Response="prop. invasive richness", site="Quarry"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.irich)$R2 %>% mutate(site="Quarry") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.irich)$Cstat  %>% mutate(Response="prop. invasive richness", site="Quarry"))

coef.summary.quarry <- coef.summary %>% unique() %>% mutate(site="Quarry")

##### STONE CANYON ######
sc.dat <- dat %>% 
  filter(site=="sc") 
rownames(sc.dat) <- 1:nrow(sc.dat) # this needs to be run in order for the standardized coefficients to be calculated

###### Soils ######
## NH4
NH4.mod.sc <- lm(NH4 ~ Burn.Severity+Treatment, data = sc.dat)
summary(NH4.mod.sc)
plot(residuals(NH4.mod.sc))
qqnorm(residuals(NH4.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NH4.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## NO3
NO3.mod.sc <- lm(NO3 ~ Burn.Severity+Treatment, data =  sc.dat)
summary(NO3.mod.sc)
plot(residuals(NO3.mod.sc))
qqnorm(residuals(NO3.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(NO3.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## PO4
P.mod.sc <- lm(PO4 ~ Burn.Severity+Treatment, data =  sc.dat)
summary(P.mod.sc)
plot(residuals(P.mod.sc))
qqnorm(residuals(P.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(P.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

###### COVER ######
## Tot
cov.mod.sc  <- lm(plot.cov ~  NO3 + NH4 +  Burn.Severity + Treatment, data = sc.dat) 
summary(cov.mod.sc)
vif(cov.mod.sc)
plot(residuals(cov.mod.sc))
qqnorm(residuals(cov.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(cov.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion invasive 
icov.mod.sc  <- lm(prop.inv.cov  ~  NO3 + NH4 + Burn.Severity + Treatment , data = sc.dat) 
summary(icov.mod.sc )
plot(residuals(icov.mod.sc ))
qqnorm(residuals(icov.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(icov.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion annual
acov.mod.sc  <- lm(prop.ann.cov  ~  NO3 + NH4  + Burn.Severity + Treatment, data = sc.dat)
summary(acov.mod.sc )
plot(residuals(acov.mod.sc ))
qqnorm(residuals(acov.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(acov.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

###### RICHNESS #####
## Total richness
rich.mod.sc <- glm.nb(plot.rich ~  NO3 + NH4 +  Burn.Severity  + Treatment, data = sc.dat) # mean(dat$plot.rich) = 20; var(dat$plot.rich) = 54 
summary(rich.mod.sc)
plot(residuals(rich.mod.sc ))
qqnorm(residuals(rich.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(rich.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion invasive richnes
irich.mod.sc <- lm(prop.inv.sp~  NO3 + NH4  + Burn.Severity + Treatment, data = sc.dat) 
summary(irich.mod.sc)
plot(residuals(irich.mod.sc ))
qqnorm(residuals(irich.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(irich.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## Proportion annual richness
arich.mod.sc <- lm(prop.ann.sp~  NO3 + NH4  + Burn.Severity  + Treatment, data = sc.dat) 
summary(arich.mod.sc)
plot(residuals(arich.mod.sc ))
qqnorm(residuals(arich.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(arich.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

###### PIECEWISE SEM BY PLANT RESPONSE #####
## total cover
modList.plot.cov <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, cov.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.cov)
coef.summary <- summary(modList.plot.cov)$coefficients
di.effect <- bind_rows(di.effect, get_effect(modList.plot.cov, target="plot.cov")$effect_long %>% mutate(Response="total cover", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.cov)$R2 %>% mutate(site="Stone Canyon"))
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.cov)$Cstat  %>% mutate(Response="total cover", site="Stone Canyon"))

## invasive cover
modList.plot.icov <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, icov.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.icov)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.icov)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.icov, target="prop.inv.cov")$effect_long %>% mutate(Response="prop.invasive cover", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.icov)$R2 %>% mutate(site="Stone Canyon")) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.icov)$Cstat  %>% mutate(Response="prop. invasive cover", site="Stone Canyon"))

## annual cover
modList.plot.acov <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, acov.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.acov)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.acov)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.acov, target="prop.ann.cov")$effect_long %>% mutate(Response="prop. annual cover", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.acov)$R2 %>% mutate(site="Stone Canyon") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.acov)$Cstat  %>% mutate(Response="prop. annual cover", site="Stone Canyon"))

## total richness
modList.plot.rich <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, rich.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.rich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.rich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.rich, target="plot.rich")$effect_long %>% mutate(Response="species richness", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.rich)$R2 %>% mutate(site="Stone Canyon") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.rich)$Cstat  %>% mutate(Response="species richness", site="Stone Canyon"))

## annual rich
modList.plot.arich <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, arich.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.arich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.arich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.arich, target="prop.ann.sp")$effect_long %>% mutate(Response="prop. annual richness", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.arich)$R2 %>% mutate(site="Stone Canyon") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.arich)$Cstat  %>% mutate(Response="prop. annual richness", site="Stone Canyon"))

## invasive rich
modList.plot.irich <- psem(NH4.mod.sc, NO3.mod.sc, P.mod.sc, irich.mod.sc, PO4 %~~% NO3, data =sc.dat)
summary(modList.plot.irich)
coef.summary <- bind_rows(coef.summary, summary(modList.plot.irich)$coefficients)
di.effect <- bind_rows(di.effect, get_effect(modList.plot.irich, target="prop.inv.sp")$effect_long %>% mutate(Response="prop. invasive richness", site="Stone Canyon"))
r2.summary <- bind_rows(r2.summary, summary(modList.plot.irich)$R2 %>% mutate(site="Stone Canyon") ) %>% unique()
Cstat.summary <- bind_rows(Cstat.summary, summary(modList.plot.irich)$Cstat  %>% mutate(Response="prop. invasive richness", site="Stone Canyon"))
coef.summary.sc <- coef.summary %>% unique() %>% mutate(site="Stone Canyon")

##### TABLES ####
coef.summary <- coef.summary.sc %>% 
  filter(!Response=="~~PO4") %>% 
  mutate(Std.Error=as.numeric(Std.Error )) %>% 
  bind_rows(., coef.summary.quarry) %>%
  set_colnames(c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value", "Std.Estimate", "Significance", "site")) %>% 
  mutate(Estimate = round_tidy(Estimate,3),
         Std.Error = round_tidy(Std.Error,3),
         Std.Estimate = round_tidy(Std.Estimate,3)) %>% 
  unite("Estimate", Estimate,  Std.Error, sep=" (") %>% 
  mutate(Estimate=paste0(Estimate, ")")) %>%
  dplyr::select(site, Response, Predictor, Estimate, Std.Estimate, Significance) %>% 
  pivot_wider(names_from=Response, values_from=c(Estimate, Std.Estimate, Significance)) %>% 
  dplyr::select(site, Predictor,
                Estimate_NH4, Std.Estimate_NH4, Significance_NH4,
                Estimate_NO3, Std.Estimate_NO3, Significance_NO3,
                Estimate_PO4, Std.Estimate_PO4, Significance_PO4,
                Estimate_plot.cov, Std.Estimate_plot.cov, Significance_plot.cov, 
                Estimate_prop.inv.cov, Std.Estimate_prop.inv.cov, Significance_prop.inv.cov,
                Estimate_prop.ann.cov, Std.Estimate_prop.ann.cov, Significance_prop.ann.cov,
                Estimate_plot.rich, Std.Estimate_plot.rich, Significance_plot.rich,
                Estimate_prop.inv.sp, Std.Estimate_prop.inv.sp, Significance_prop.inv.sp,
                Estimate_prop.ann.sp, Std.Estimate_prop.ann.sp, Significance_prop.ann.sp) 
coef.summary %>% 
  arrange(site) %>% 
  write.csv(here("output", "SEM-coefsummary.csv"), row.names=F)

r2.summary %>% 
  write.csv(here("output", "SEM-r2summary.csv"), row.names=F)

##### Figure 3 #####
p1 <- estimate_relation(NH4.mod.q, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) +
    geom_line()+geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
    geom_point(data=q.dat, aes(x=Burn.Severity, y=NH4, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) +
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) +
    xlab("Burn Severity") +
    ylab(expression(paste("soil NH"[4], " (mg/kg)"))) + 
    theme(legend.position = "bottom", legend.title=element_blank())+ylim(0,max(dat$NH4)+0.05)


p2 <- estimate_relation(NO3.mod.q, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) + 
    geom_line() + 
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=q.dat, aes(x=Burn.Severity, y=NO3, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) + 
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) + 
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) +
    xlab("Burn Severity") + 
    ylab(expression(paste("soil NO"[3], " (mg/kg)"))) +
    theme(legend.position = "bottom", legend.title=element_blank())

p3 <- estimate_relation(P.mod.q, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) + 
    geom_line() + 
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=q.dat, aes(x=Burn.Severity, y=PO4, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) + 
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    xlab("Burn Severity") + 
    ylab(expression(paste("soil PO"[4], " (mg/kg)"))) + 
    theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "bottom", legend.title=element_blank())

p4 <- estimate_relation(NH4.mod.sc, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) + 
    geom_line() + 
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
    geom_point(data=sc.dat, aes(x=Burn.Severity, y=NH4, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) + 
    scale_x_continuous(breaks = c(1, 2), labels = c("unburned", "low")) + 
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    xlab("Burn Severity") +
    ylab(expression(paste("soil NH"[4], " (mg/kg)"))) + 
    theme(legend.position = "bottom", legend.title=element_blank())+ylim(0,max(dat$NH4)+0.05)


p5 <- estimate_relation(NO3.mod.sc, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) + 
    geom_line() + 
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
    geom_point(data=sc.dat, aes(x=Burn.Severity, y=NO3, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) + 
    scale_x_continuous(breaks = c(1, 2), labels = c("unburned", "low")) + 
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    xlab("Burn Severity") +
    ylab(expression(paste("soil NO"[3], " (mg/kg)"))) +
    theme(plot.title = element_text(hjust = 0.5),legend.position = "bottom", legend.title=element_blank())

p6 <- estimate_relation(P.mod.sc, by=c("Burn.Severity", "Treatment")) %>% 
  as.data.frame() %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>% 
  mutate(Predicted=ifelse(Predicted<0, 0, Predicted)) %>% 
  ggplot(aes(x=Burn.Severity, y=Predicted, col=factor(Treatment), fill=factor(Treatment))) + 
    geom_line() + 
    geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
    scale_x_continuous(breaks = c(1, 2), labels = c("unburned", "low")) + 
    scale_color_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant")) + 
    scale_fill_manual(values=c("gray25","darkred"), labels=c("control", "Fire retardant") )+ 
    xlab("Burn Severity") + 
    ylab(expression(paste("soil PO"[4], " (mg/kg)"))) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    geom_point(data=sc.dat, aes(x=Burn.Severity, y=PO4, col=factor(Treatment)), position=position_jitter(width=0.1), alpha=0.75) + 
    theme(legend.position = "bottom", legend.title=element_blank())

p123456 <- p1 +p2 +p3+p4 +p5 +p6+plot_layout(guides="collect",  nrow=3, byrow=F, widths=c(2,1), axis_titles="collect")+ plot_annotation(tag_levels = 'a') & theme(legend.position = "bottom", legend.title=element_blank())

ggsave(p123456, filename=here("output", "fig3-soils.jpg"), width =3.5, height = 6, units = "in", dpi = 300)




##### FIGURE 4 ####
p1 <- bind_rows(estimate_relation(cov.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(cov.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.cov),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("total cover (%)") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p2 <- bind_rows(estimate_relation(cov.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(cov.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.cov),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("total cover (%)") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p3 <- estimate_relation(cov.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry") %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.cov), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("total cover (%)") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p4 <- bind_rows(estimate_relation(cov.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(cov.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("total cover (%)") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 
  
p5 <- bind_rows(estimate_relation(cov.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(cov.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("total cover (%)") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

### proportion non-native
p6 <- bind_rows(estimate_relation(icov.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(icov.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.cov),alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\ncover") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p7 <- bind_rows(estimate_relation(icov.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(icov.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.cov), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\ncover") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p8 <- estimate_relation(icov.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry") %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.cov), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\ncover") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p9 <- bind_rows(estimate_relation(icov.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(icov.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\ncover") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 

p10<- bind_rows(estimate_relation(icov.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(icov.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\ncover") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

### proportion annual 
p11 <- bind_rows(estimate_relation(acov.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(acov.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.cov), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\ncover") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p12 <- bind_rows(estimate_relation(acov.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(acov.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.cov),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\ncover") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p13 <- estimate_relation(acov.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry") %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.cov),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\ncover") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p14 <- bind_rows(estimate_relation(acov.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(acov.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\ncover") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 

p15 <- bind_rows(estimate_relation(acov.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(acov.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc"))%>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.cov), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\ncover") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

pcover <- p1 + p2 +p3 +p4 +p5 +p6 +p7 +p8 +p9 +p10 +p11 +p12 +p13 +p14 +p15 +plot_layout(guides="collect", nrow=5, byrow=F)+ plot_annotation(tag_levels = 'a') & theme(legend.position = "bottom", legend.title=element_blank()) 
pcover %>% ggsave(filename=here("output", "fig4-cover.jpg"), width = 7, height = 8, units = "in", dpi = 300)

##### FIGURE 5 ####
p1 <- bind_rows(estimate_relation(rich.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(rich.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.rich),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("species richness") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p2 <- bind_rows(estimate_relation(rich.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(rich.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.rich),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("species richness") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p3 <- estimate_relation(rich.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry") %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.rich), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("species richness") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p4 <- bind_rows(estimate_relation(rich.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(rich.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.rich), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("species richness") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 

p5 <- bind_rows(estimate_relation(rich.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(rich.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=plot.rich), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("species richness") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

### proportion non-native
p6 <- bind_rows(estimate_relation(irich.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(icov.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.sp),alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\nspecies") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p7 <- bind_rows(estimate_relation(irich.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(irich.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.sp), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\nspecies") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p8 <- estimate_relation(irich.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry") %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.sp), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\nspecies") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p9 <- bind_rows(estimate_relation(irich.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(irich.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.sp), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\nspecies") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 

p10<- bind_rows(estimate_relation(irich.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(irich.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.inv.sp), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. non-native\nspecies") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

### proportion annual 
p11 <- bind_rows(estimate_relation(arich.mod.q, by="NH4") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(arich.mod.sc, by="NH4")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NH4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.sp), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\nspecies") + xlab(expression(paste("soil NH"[4], " (mg/g)"))) 

p12 <- bind_rows(estimate_relation(arich.mod.q, by="NO3") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(arich.mod.sc, by="NO3")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=NO3, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.sp),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\nspecies") + xlab(expression(paste("soil NO"[3], " (mg/g)"))) 

p13 <- estimate_relation(arich.mod.q, by="PO4") %>% as.data.frame() %>% mutate(site="quarry")%>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=PO4, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.sp),  alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\nspecies") + xlab(expression(paste("soil PO"[4], " (mg/g)"))) 

p14 <- bind_rows(estimate_relation(arich.mod.q, by="Burn.Severity") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(arich.mod.sc, by="Burn.Severity")  %>% as.data.frame()  %>% mutate(site="sc")) %>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Burn.Severity, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.sp), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\nspecies") + xlab("burn severity") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) 

p15 <- bind_rows(estimate_relation(arich.mod.q, by="Treatment") %>% as.data.frame() %>% mutate(site="quarry"), estimate_relation(arich.mod.sc, by="Treatment")  %>% as.data.frame()  %>% mutate(site="sc"))%>% 
  mutate(CI_low=ifelse(CI_low<0, 0, CI_low)) %>%
  ggplot(aes(x=Treatment, y=Predicted, col=site, fill=site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.25, color=NA) + 
  geom_point(data=dat, aes(y=prop.ann.sp), position=position_jitter(width=0.1), alpha=0.75) + 
  scale_color_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  scale_fill_manual(values=c("#018571", "#A6611A"), labels=c("Quarry", "Stone Canyon")) +
  ylab("prop. annual\nspecies") + xlab("fire retardant") +
  scale_x_continuous(breaks = c(0, 1), labels = c("control", "treated")) 

prichness <- p1 + p2 +p3 +p4 +p5 +p6 +p7 +p8 +p9 +p10 +p11 +p12 +p13 +p14 +p15 +plot_layout(guides="collect", nrow=5, byrow=F)+ plot_annotation(tag_levels = 'a') & theme(legend.position = "bottom", legend.title=element_blank()) 
prichness %>% ggsave(filename=here("output", "fig5-richness.jpg"), width = 7, height = 8, units = "in", dpi = 300)

##### SUPPLEMENTAL FIGURES #####
###### FIGURE S2 - DIRECT AND INDIRECT EFFECTS #####
ps2 <- di.effect %>% 
  mutate(Variable=factor(Variable, levels=c("Burn.Severity", "Treatment", "Burn.Severity:Treatment", "NH4", "NO3", "PO4" ), labels=c("Burn severity", "Fire retardant", "Burn Severity:fire retardant", "soil NH4", "soil NO3", "soil PO4")),
         Effect_Type=factor(Effect_Type, levels=c("Direct_Effect", "Indirect_Effect", "Total_Effect"), labels=c("Direct", "Indirect", "Total")),
         Response=factor(Response, levels=c("total cover", "prop.invasive cover", "prop. annual cover", "species richness",  "prop. invasive richness","prop. annual richness"), labels=c("total cover", "prop. non-native cover", "prop. annual cover", "species richness", "prop. non-native species", "prop. annual species"), ordered=T)) %>%
  ggplot(aes(x=Variable, y=Effect_Value, fill=Effect_Type)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(Response~site, , labeller = label_wrap_gen(width = 16))  + 
  xlab("") + ylab("standardized effect size") + 
  theme(legend.position = "bottom", legend.title=element_blank())+coord_flip()

ggsave(ps2, filename=here("output", "figS2-direct-indirect-effects.jpg"), width = 7, height = 7, units = "in", dpi = 300)



###### FIGURE S5 - CORRELATION PLOT - QUARRY #####
library(ggcorrplot)
p1 <- q.dat[,c("NH4", "NO3", "PO4", "plot.cov", "prop.inv.cov", "prop.ann.cov", "plot.rich", "prop.inv.sp", "prop.ann.sp")] %>% 
  set_colnames(c("NH4", "NO3", "PO4", "total cover", "prop. non-native cover", "prop. annual cover", "species richness", "prop. non-native species", "prop. annual species")) %>% 
    cor( use="pairwise.complete.obs") %>%
  ggcorrplot(method="square", type="lower", lab=T, lab_size=3) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(title="Quarry")

ggsave(p1, filename=here("output", "figS5-correlation-matrix-Quarry.jpg"), width = 5, height = 5, units = "in", dpi = 300)


###### FIGURE S5 - CORRELATION PLOT - STONE CANYON #####
p2 <- sc.dat[,c("NH4", "NO3", "PO4", "plot.cov", "prop.inv.cov", "prop.ann.cov", "plot.rich", "prop.inv.sp", "prop.ann.sp")] %>% 
  set_colnames(c("NH4", "NO3", "PO4", "total cover", "prop. non-native cover", "prop. annual cover", "species richness", "prop. non-native species", "prop. annual species")) %>% 
  cor( use="pairwise.complete.obs") %>%
  ggcorrplot(method="square", type="lower", lab=T, lab_size=3) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(title="Stone Canyon")
ggsave(p2, filename=here("output", "figS6-correlation-matrix-StoneCanyon.jpg"), width = 5, height = 5, units = "in", dpi = 300)

###### FIGURE S6 - CONTRIBUTION OF NATIVES TO ANNUAL COVER #####
dat.cover.annuals.fr <- bind_rows(comm.dat.full.q, comm.dat.full.sc) %>% 
  left_join(., env[,c("plot", "site", "trt", "sev")], by = c("plot" = "plot")) %>% 
  left_join(., sp.info, by=c("code"="code")) %>% 
  filter(duration=="annual") %>%
  group_by(plot, site, trt, sev, duration, status) %>% 
  summarize(cov=sum(plotavg)) %>% 
  pivot_wider(names_from=status, values_from=cov) %>% 
  mutate(p.annual=N/(I+N)) %>% 
  group_by(site, trt, sev) %>% 
  mutate(sev=as.numeric(factor(sev, levels=c("unburn", "low", "mod", "high"), labels=1:4, ordered=T)),
         site=factor(site, levels=c("quarry", "sc"), labels=c("Quarry", "Stone Canyon"))) %>%
  ggplot(aes(x=sev, y=p.annual, col=trt))+
  geom_point(stat = "summary", fun = "mean", position=position_dodge(width=0.5), shape=17, size=3, geom='pointrange') +  # Plots the mean point
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5))+facet_wrap(.~site, nrow=1, scales="free_x", space="free_x") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) +
  scale_color_manual(values=c("gray25","red"), labels=c("control", "Fire retardant"))+theme(legend.position = "bottom", legend.title=element_blank())+xlab("burn severity")+ylab("contribution of natives to annual cover")
ggsave(dat.cover.annuals.fr, filename=here("output", "figS6-prop-annual-cover-by-severity-treatment.jpg"), width = 5, height = 5, units = "in", dpi = 300)






##### SUPPLEMENTAL TABLES #####
###### TABLE S1 - COVER AND RICHNESS IN UNBURNED AND UNTREATED PLOTS #####
S1 <- dat %>% 
  filter(Treatment==0 & Burn.Severity==1) %>% 
  group_by(site) %>% 
  summarise(native = mean(plot.cov-invasive.cov, na.rm=T), 
            native.sd = sd(plot.cov-invasive.cov, na.rm=T),
            invasive = mean(invasive.cov, na.rm=T), 
            invasive.sd = sd(invasive.cov, na.rm=T), 
            annual=mean(annual.cov), 
            annual.sd=sd(annual.cov), 
            native.spp = mean(plot.rich-invasive.rich), 
            native.spp.sd =sd(plot.rich-invasive.rich), 
            invasive.spp = mean(invasive.rich), 
            invasive.spp.sd = sd(invasive.rich), 
            annual.spp = mean(annual.rich),
            annual.spp.sd = sd(annual.rich),
            ) %>% 
  mutate(across(where(is.numeric), round_tidy, 1)) %>%
  unite("native", native, native.sd, sep=" ± ") %>% 
  unite("invasive", invasive, invasive.sd, sep=" ± ") %>%
  unite("annual", annual, annual.sd, sep=" ± ") %>%
  unite("native.spp", native.spp, native.spp.sd, sep=" ± ") %>% 
  unite("invasive.spp", invasive.spp, invasive.spp.sd, sep=" ± ") %>%
  unite("annual.spp", annual.spp, annual.spp.sd, sep=" ± ") 
S1 %>% 
  write_csv(here("output", "S1-Veg-UntreatedUnburned.csv"))

###### TABLE S2 - COMMON SPECIES #####
n.plots <- env %>% group_by(site, trt, sev) %>% summarise(n.plots=n(), .groups = "drop")

com.dat <- comm.sum.long %>% 
  left_join(., env[,c("plot", "site")], by = c("plot" = "plot")) %>% 
  mutate(presence=1)

exq <- expand.grid(com.dat %>% filter(site=="quarry") %>% pull(plot) %>% unique(), com.dat %>% filter(site=="quarry") %>% pull(code) %>% unique()) %>% 
  set_colnames(c("plot", "code")) 

comm.dat.full.q <- com.dat[com.dat$site=="quarry", c("plot", "code", "plotavg", "presence")] %>% 
  right_join(., exq,  by=c("plot"="plot", "code"="code")) %>% 
  mutate(plotavg=ifelse(is.na(plotavg), 0, plotavg),
    presence=ifelse(is.na(presence), 0, presence)) 

exsc <- expand.grid(com.dat %>% filter(site=="sc") %>% pull(plot) %>% unique(), com.dat %>% filter(site=="sc") %>% pull(code) %>% unique()) %>% 
  set_colnames(c("plot", "code")) 

comm.dat.full.sc <- com.dat[com.dat$site=="sc",c("plot", "code", "plotavg", "presence")] %>% 
  right_join(., exsc,  by=c("plot"="plot", "code"="code")) %>% 
  mutate(plotavg=ifelse(is.na(plotavg), 0, plotavg),
         presence=ifelse(is.na(presence), 0, presence)) 

dat.cover <- bind_rows(comm.dat.full.q, comm.dat.full.sc) %>% 
  left_join(., env[,c("plot", "site", "trt", "sev")], by = c("plot" = "plot"))  %>% 
  filter(trt=="con" & sev=="unburn") %>% 
  group_by(site, code) %>% 
  summarise(mean = mean(plotavg, na.rm=T), prop=sum(presence)/n(), .groups = "drop")  %>% 
  left_join(., sp.info, by=c("code"="code"))

s2 <- dat.cover %>% 
  group_by(site) %>% 
 arrange(desc(prop)) %>% slice(1:10) %>% 
  dplyr::select(site, species, prop,  mean, duration, functional.group, status) %>% 
  mutate(prop=round_tidy(prop, 2), mean=round_tidy(mean, 1))
s2 %>% write_csv(here("output", "S2-CommonSpp-UntreatedUnburned.csv"))


### common species in burned
highsev.spp <- bind_rows(comm.dat.full.q, comm.dat.full.sc) %>% 
  left_join(., env[,c("plot", "site", "trt", "sev")], by = c("plot" = "plot"))  %>% 
  filter(trt=="con" & sev=="high") %>% 
  group_by(site, code) %>% 
  summarise(mean = mean(plotavg, na.rm=T), prop=sum(presence)/n(), .groups = "drop")  %>% 
  left_join(., sp.info, by=c("code"="code")) %>% 
  group_by(site) %>% 
  arrange(desc(mean)) %>% slice(1:10) %>% 
  dplyr::select(site, species, prop,  mean, duration, functional.group, status) %>% 
  mutate(prop=round_tidy(prop, 2), mean=round_tidy(mean, 1))


### common species in treated
highsev.spp <- bind_rows(comm.dat.full.q, comm.dat.full.sc) %>% 
  left_join(., env[,c("plot", "site", "trt", "sev")], by = c("plot" = "plot"))  %>% 
  filter(trt=="fr") %>% 
  group_by(site, code) %>% 
  summarise(mean = mean(plotavg, na.rm=T), prop=sum(presence)/n(), .groups = "drop")  %>% 
  left_join(., sp.info, by=c("code"="code")) %>% 
  group_by(site) %>% 
  arrange(desc(mean)) %>% slice(1:10) %>% 
  dplyr::select(site, species, prop,  mean, duration, functional.group, status) %>% 
  mutate(prop=round_tidy(prop, 2), mean=round_tidy(mean, 1))


###### FIGURE S6 - CONTRIBUTION OF NON-NATIVE ANNUALS TO TOTAL COVER #####
dat.cover.annuals.fr <- bind_rows(comm.dat.full.q, comm.dat.full.sc) %>% 
  left_join(., env[,c("plot", "site", "trt", "sev")], by = c("plot" = "plot")) %>% 
  left_join(., sp.info, by=c("code"="code")) %>% 
  filter(duration=="annual") %>%
  group_by(plot, site, trt, sev, duration, status) %>% 
  summarize(cov=sum(plotavg)) %>% 
  pivot_wider(names_from=status, values_from=cov) %>%
  mutate(sev=as.numeric(factor(sev, levels=c("unburn", "low", "mod", "high"), labels=1:4, ordered=T)),
         site=factor(site, levels=c("quarry", "sc"), labels=c("Quarry", "Stone Canyon")),
         p.native = N/(I+N)) %>%
  ggplot(aes(x=sev, y=p.native, col=trt))+
  geom_point(stat = "summary", fun = "mean", position=position_dodge(width=0.5), shape=17, size=3, geom='pointrange') +  # Plots the mean point
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5))+facet_wrap(.~site, nrow=1, scales="free_x", space="free_x") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("unburned", "low", "moderate", "high")) +
  scale_color_manual(values=c("gray25","red"), labels=c("control", "Fire retardant"))+theme(legend.position = "bottom", legend.title=element_blank())+xlab("burn severity")+ylab("contribution of natives to annual cover")
ggsave(dat.cover.annuals.fr, filename=here("output", "figS6-prop-annual-cover-by-severity-treatment.jpg"), width = 5, height = 3, units = "in", dpi = 300)




###### TABLE S4 - OBSERVED NUTRIENT VALUES BY BURN SEVERITY AND TREATMENT #####
S4a <- dat %>% 
  group_by(site, Burn.Severity, Treatment) %>%
  summarise(NH4 = mean(NH4, na.rm=T), 
            NO3 = mean(NO3, na.rm=T), 
            PO4 = mean(PO4, na.rm=T), 
            pH = mean(pH, na.rm=T), 
            Na = mean(Na, na.rm=T),
            K = mean(K, na.rm=T),
            SO4 = mean(SO4, na.rm=T),
            .groups = "drop") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=c(1, 2, 3, 4), labels=c("unburned", "low", "moderate", "high"), ordered=T),Treatment=factor(Treatment, levels=c(0, 1), labels=c("control", "treated"))) %>%
  pivot_longer(cols=c(NH4, NO3, PO4, pH,  Na, K,  SO4), names_to="nutrient", values_to="value") %>% 
  pivot_wider(names_from = Treatment, values_from = value) %>% 
  arrange(site,  nutrient, Burn.Severity) %>% 
  mutate(adifference=round_tidy(treated-control, 3), 
         pchange = round_tidy(((treated-control)/control*100), 0)) %>%  
  dplyr::select(-control, -treated)

S4b <- dat %>% 
  group_by(site, Burn.Severity, Treatment) %>%
  summarise(NH4.mean=mean(NH4, na.rm=T), 
            NH4.sd=sd(NH4, na.rm=T), 
            NO3.mean=mean(NO3, na.rm=T), 
            NO3.sd=sd(NO3, na.rm=T), 
            PO4.mean=mean(PO4, na.rm=T), 
            PO4.sd=sd(PO4, na.rm=T),
            pH.mean = mean(pH, na.rm=T), 
            pH.sd = sd(pH, na.rm=T),
            Na.mean = mean(Na, na.rm=T),
            Na.sd = sd(Na, na.rm=T),
            K.mean = mean(K, na.rm=T),
            K.sd = sd(K, na.rm=T),
            SO4.mean = mean(SO4, na.rm=T),
            SO4.sd = sd(SO4, na.rm=T), 
            .groups = "drop") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=c(1, 2, 3, 4), labels=c("unburned", "low", "moderate", "high"), ordered=T),
         Treatment=factor(Treatment, levels=c(0, 1), labels=c("control", "treated")),
         pH.mean=round_tidy(pH.mean, 2),
         pH.sd=round_tidy(pH.sd, 2)) %>% 
  mutate(across(where(is.numeric), round_tidy, 3)) %>%
  unite("NH4", NH4.mean, NH4.sd, sep=" ± ") %>%
  unite("NO3", NO3.mean, NO3.sd, sep=" ± ") %>%
  unite("PO4", PO4.mean, PO4.sd, sep=" ± ") %>%
  unite("pH", pH.mean, pH.sd, sep=" ± ")  %>%
  unite("Na", Na.mean, Na.sd, sep=" ± ") %>%
  unite("K", K.mean, K.sd, sep=" ± ") %>%
  unite("SO4", SO4.mean, SO4.sd, sep=" ± ") %>%
  pivot_longer(cols=c(NH4, NO3, PO4, pH, Na, K, SO4),names_to="nutrient", values_to="value") %>% 
  pivot_wider(names_from = Treatment, values_from = value) %>% 
  arrange(site,  nutrient, Burn.Severity)

S4ab <- full_join(S4b, S4a, by=c("site", "Burn.Severity", "nutrient")) %>% 
  arrange(site, nutrient, Burn.Severity) 

S4ab %>% 
  write_csv(here("output", "S4-observed-soilsXTreatmentXBurn.csv"))


###### TABLE S5 - PARAMAETER ESTIMATES FOR OTHER SOIL NUTRIENTS AND CHEMICALS #####
# QUARRY
## pH
pH.mod.q <- lm(pH ~ Burn.Severity*Treatment, data = q.dat)
summary(pH.mod.q)
plot(residuals(pH.mod.q))
qqnorm(residuals(pH.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(pH.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## Na
Na.mod.q <- lm(Na ~ Burn.Severity*Treatment, data =  q.dat)
summary(Na.mod.q)
plot(residuals(Na.mod.q))
qqnorm(residuals(Na.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(Na.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## K
K.mod.q <- lm(K ~ Burn.Severity*Treatment, data =  q.dat)
summary(K.mod.q)
plot(residuals(K.mod.q))
qqnorm(residuals(K.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(K.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

## SO4
SO4.mod.q <- lm(SO4 ~ Burn.Severity*Treatment, data =  q.dat)
summary(SO4.mod.q)
plot(residuals(SO4.mod.q))
qqnorm(residuals(SO4.mod.q), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(SO4.mod.q), col = "red", lwd = 2, lty = 2) ## okay enough

modList.other <- psem(pH.mod.q, Na.mod.q, K.mod.q, SO4.mod.q, data =q.dat)
coef.summary.other <- summary(modList.other )$coefficients %>% 
  as.data.frame() %>%  
  set_colnames(c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.value", "Std.Estimate", "Significance")) %>%
  mutate(site="Quarry") 
r2.summary.other  <- summary(modList.other )$R2 %>% as.data.frame() %>% mutate(site="Quarry")


# STONE CANYON
## pH
pH.mod.sc <- lm(pH ~ Burn.Severity+Treatment, data = sc.dat)
summary(pH.mod.sc)
plot(residuals(pH.mod.sc))
qqnorm(residuals(pH.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(pH.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## Na
Na.mod.sc <- lm(Na ~ Burn.Severity+Treatment, data =  sc.dat)
summary(Na.mod.sc)
plot(residuals(Na.mod.sc))
qqnorm(residuals(Na.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(Na.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## K
K.mod.sc <- lm(K ~ Burn.Severity+Treatment, data =  sc.dat)
summary(K.mod.sc)
plot(residuals(K.mod.sc))
qqnorm(residuals(K.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(K.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

## SO4
SO4.mod.sc <- lm(SO4 ~ Burn.Severity+Treatment, data =  sc.dat)
summary(SO4.mod.sc)
plot(residuals(SO4.mod.sc))
qqnorm(residuals(SO4.mod.sc), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(SO4.mod.sc), col = "red", lwd = 2, lty = 2) ## okay enough

modList.other <- psem(pH.mod.sc, Na.mod.sc, K.mod.sc, SO4.mod.sc, data =sc.dat)
coef.summary.other <- bind_rows(coef.summary.other, 
                                summary(modList.other )$coefficients %>% 
                                  as.data.frame() %>%  
                                  set_colnames(c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.value", "Std.Estimate", "Significance")) %>%
                                  mutate(site="Stone Canyon"))
r2.summary.other  <- bind_rows(r2.summary.other, summary(modList.other )$R2 %>% as.data.frame() %>% mutate(site="Stone Canyon"))
r2.summary.other %>% write.csv(here("output", "SEM-r2summary-other.csv"), row.names=F)

coef.summary.other  %>% 
  mutate(Estimate=round_tidy(Estimate, 3), 
         Std.Error=round_tidy(Std.Error, 3), 
         Crit.Value=round_tidy(Crit.Value, 3), 
         P.value=round_tidy(P.value, 3), 
         Std.Estimate=round_tidy(Std.Estimate, 3)) %>% 
           unite("Estimate", Estimate,  Std.Error, sep=" (") %>% 
  mutate(Estimate=paste0(Estimate, ")")) %>%
  dplyr::select(site, Response, Predictor, Estimate, Std.Estimate, Significance) %>% 
  pivot_wider(names_from=Response, values_from=c(Estimate, Std.Estimate, Significance)) %>% 
  dplyr::select(site, Predictor,
                Estimate_pH, Std.Estimate_pH, Significance_pH,
                Estimate_Na, Std.Estimate_Na, Significance_Na,
                Estimate_K, Std.Estimate_K, Significance_K,
                Estimate_SO4, Std.Estimate_SO4, Significance_SO4) %>% 
   write.csv(here("output", "SEM-coefsummary-other.csv"), row.names=F)


###### TABLE S6 - MARGINAL MEANS OF SOIL NUTRIENTS BY TREATMENT AND FIRE SEVERITY #####
S6.q <-  list (NH4.mod.q, NO3.mod.q, P.mod.q, pH.mod.q, Na.mod.q, K.mod.q, SO4.mod.q) %>% 
  set_names(c("NH4", "NO3", "PO4", "pH", "Na", "K", "SO4")) %>% 
  lapply(estimate_means, by=c("Treatment", "Burn.Severity")) %>% 
  bind_rows(.id="response") %>% 
  mutate(Mean=ifelse(Mean<0, 0, Mean)) %>% 
  dplyr::select(response, Treatment,Burn.Severity, Mean, SE) %>%
  pivot_wider(names_from = c(Treatment), values_from = c(Mean, SE)) %>% 
  mutate(adif=round_tidy(`Mean_1`-`Mean_0`, 3),
         pchange= round_tidy((`Mean_1`-`Mean_0`)/`Mean_0`*100,0)) %>% 
  mutate(across(where(is.numeric), round_tidy, 3)) %>% 
  unite(Control, Mean_0, SE_0, sep=" ± ") %>% 
  unite(Treated, Mean_1, SE_1, sep=" ± ") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=round_tidy(1:4, 3), labels=c("unburned", "low", "moderate","high")), 
         pchange=ifelse( pchange==Inf, NA,  pchange), 
         site="Quarry")

S6.sc <-  list (NH4.mod.sc, NO3.mod.sc, P.mod.sc, pH.mod.sc, Na.mod.sc, K.mod.sc, SO4.mod.sc) %>% 
  set_names(c("NH4", "NO3", "PO4", "pH", "Na", "K", "SO4")) %>% 
  lapply(estimate_means, by=c("Treatment", "Burn.Severity")) %>% 
  bind_rows(.id="response") %>% 
  mutate(Mean=ifelse(Mean<0, 0, Mean)) %>% 
  dplyr::select(response, Treatment,Burn.Severity, Mean, SE) %>%
  pivot_wider(names_from = c(Treatment), values_from = c(Mean, SE)) %>% 
  mutate(adif=round_tidy(`Mean_1`-`Mean_0`, 3),
         pchange= round_tidy((`Mean_1`-`Mean_0`)/`Mean_0`*100,0)) %>% 
  mutate(across(where(is.numeric), round_tidy, 3)) %>% 
  unite(Control, Mean_0, SE_0, sep=" ± ") %>% 
  unite(Treated, Mean_1, SE_1, sep=" ± ") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=round_tidy(1:2, 3), labels=c("unburned", "low")), 
         pchange=ifelse( pchange==Inf, NA,  pchange), 
         site="Stone Canyon")

S6 <- bind_rows(S6.q, S6.sc) %>% 
  write.csv(here("output", "S6-marginalmeans-soilsXTreatmentXBurn.csv"), row.names=F)



###### TABLE S7 - PLANT COVER BY BURN SEVERITY AND TREATMENT #####
S7a <- dat %>% 
  group_by(site, Burn.Severity, Treatment) %>%
  summarise(plot.cov = mean(plot.cov, na.rm=T), 
            prop.inv.cov = mean(prop.inv.cov , na.rm=T), 
            prop.ann.cov  = mean(prop.ann.cov, na.rm=T), 
            plot.rich = mean(plot.rich, na.rm=T), 
            prop.inv.sp = mean(prop.inv.sp, na.rm=T), 
            prop.ann.sp  = mean(prop.ann.sp, na.rm=T), 
            .groups = "drop") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=c(1, 2, 3, 4), labels=c("unburned", "low", "moderate", "high"), ordered=T),Treatment=factor(Treatment, levels=c(0, 1), labels=c("control", "treated"))) %>%
  pivot_longer(cols=c(plot.cov, prop.inv.cov, prop.ann.cov, plot.rich, prop.inv.sp, prop.ann.sp), names_to="community", values_to="value") %>% 
  pivot_wider(names_from = Treatment, values_from = value) %>% 
  arrange(site,  community, Burn.Severity) %>% 
  mutate(difference=round_tidy(treated-control,2), 
         pchange =round(((treated-control)/control)*100, 0)) %>% 
  dplyr::select(-control, -treated)

S7b <- dat %>% 
  group_by(site, Burn.Severity, Treatment) %>%
  summarise(plot.cov.mean = mean(plot.cov, na.rm=T), 
            plot.cov.sd =sd(plot.cov, na.rm=T),
            prop.inv.cov.mean  = mean(prop.inv.cov , na.rm=T), 
            prop.inv.cov.sd = sd(prop.inv.cov , na.rm=T),
            prop.ann.cov.mean   = mean(prop.ann.cov, na.rm=T), 
            prop.ann.cov.sd = sd(prop.ann.cov, na.rm=T),
            plot.rich.mean  = mean(plot.rich, na.rm=T), 
            plot.rich.sd = sd(plot.rich, na.rm=T),
            prop.inv.sp.mean  = mean(prop.inv.sp, na.rm=T), 
            prop.inv.sp.sd = sd(prop.inv.sp, na.rm=T),
            prop.ann.sp.mean   = mean(prop.ann.sp, na.rm=T), 
            prop.ann.sp.sd = sd(prop.ann.sp, na.rm=T),
            .groups = "drop")  %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=c(1, 2, 3, 4), labels=c("unburned", "low", "moderate", "high"), ordered=T),
         Treatment=factor(Treatment, levels=c(0, 1), labels=c("control", "treated"))) %>% 
  mutate(plot.cov.mean=round_tidy(plot.cov.mean, 1), 
         plot.cov.sd  = round_tidy(plot.cov.sd, 1),
         plot.rich.mean=round_tidy(plot.rich.mean, 1),
         plot.rich.sd =round_tidy(plot.rich.sd, 1)) %>% 
  mutate(across(where(is.numeric), round_tidy, 2)) %>%
  unite("plot.cov", plot.cov.mean ,plot.cov.sd, sep=" ± ") %>%
  unite("prop.inv.cov", prop.inv.cov.mean , prop.inv.cov.sd, sep=" ± ") %>%
  unite("prop.ann.cov", prop.ann.cov.mean , prop.ann.cov.sd, sep=" ± ") %>%
  unite("plot.rich", plot.rich.mean , plot.rich.sd, sep=" ± ") %>%
  unite("prop.inv.sp", prop.inv.sp.mean , prop.inv.sp.sd, sep=" ± ") %>%
  unite("prop.ann.sp", prop.ann.sp.mean , prop.ann.sp.sd, sep=" ± ") %>%
  pivot_longer(cols=c(plot.cov, prop.inv.cov, prop.ann.cov, plot.rich, prop.inv.sp, prop.ann.sp), names_to="community", values_to="value") %>% 
  pivot_wider(names_from = Treatment, values_from = value) %>% 
  arrange(site,  community, Burn.Severity)

S7ab <- full_join(S7b, S7a, by=c("site", "Burn.Severity", "community")) %>% 
  arrange(site, community, Burn.Severity) 

S7ab %>% 
  write_csv(here("output", "S7-obseved-vegXTreatmentXBurn.csv"))


###### TABLE S8 - MARGINAL MEAN PLANT COVER BY TREATMENT AND BURN SEVERITY #####

S8.q <- list (cov.mod.q, icov.mod.q, acov.mod.q, rich.mod.q, irich.mod.q, arich.mod.q) %>% 
  set_names(c("total cover", "prop. non-native cover", "prop. annual cover", "species richnes", "prop. non-native richness", "prop. annual richness")) %>% 
  lapply(estimate_means, by=c("Treatment", "Burn.Severity")) %>% 
  bind_rows(.id="response") %>% 
  mutate(Mean=ifelse(Mean<0, 0, Mean)) %>% 
  dplyr::select(response, Treatment,Burn.Severity, Mean, SE) %>%
  pivot_wider(names_from = c(Treatment), values_from = c(Mean, SE)) %>% 
  mutate(adif=round_tidy(`Mean_1`-`Mean_0`, 2),
         pchange= round_tidy((`Mean_1`-`Mean_0`)/`Mean_0`*100,0)) %>% 
  mutate(across(where(is.numeric), round_tidy, 2)) %>% 
  unite(Control, Mean_0, SE_0, sep=" ± ") %>% 
  unite(Treated, Mean_1, SE_1, sep=" ± ") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=round_tidy(1:4, 2), labels=c("unburned", "low", "moderate","high")), 
         pchange=ifelse( pchange==Inf, NA,  pchange), 
         site="Quarry")

S8.sc <- list (cov.mod.sc, icov.mod.sc, acov.mod.sc, rich.mod.sc, irich.mod.sc, arich.mod.sc) %>% 
  set_names(c("total cover", "prop. non-native cover", "prop. annual cover", "species richnes", "prop. non-native richness", "prop. annual richness")) %>% 
  lapply(estimate_means, by=c("Treatment", "Burn.Severity")) %>% 
  bind_rows(.id="response") %>% 
  mutate(Mean=ifelse(Mean<0, 0, Mean)) %>% 
  dplyr::select(response, Treatment,Burn.Severity, Mean, SE) %>%
  pivot_wider(names_from = c(Treatment), values_from = c(Mean, SE)) %>% 
  mutate(adif=round_tidy(`Mean_1`-`Mean_0`, 2),
         pchange= round_tidy((`Mean_1`-`Mean_0`)/`Mean_0`*100,0)) %>% 
  mutate(across(where(is.numeric), round_tidy, 2)) %>% 
  unite(Control, Mean_0, SE_0, sep=" ± ") %>% 
  unite(Treated, Mean_1, SE_1, sep=" ± ") %>% 
  mutate(Burn.Severity= factor(Burn.Severity, levels=round_tidy(1:4, 2), labels=c("unburned", "low", "moderate","high")), 
         pchange=ifelse( pchange==Inf, NA,  pchange), 
         site="Stone Canyon")

S8 <- bind_rows(S8.q, S8.sc) %>% 
  write.csv(here("output", "S8-marginalmeans-vegXTreatmentXBurn.csv"), row.names=F)



###### SUPPLEMENTAL ORDINATIONS #####

pal <- c("#f6d746", "#e55c30", "#84206b", "#140b34")
pal2 <- c("black", "red")


## Quarry Fire  
Quarry <- as.matrix(comm.sum[match(env$plot[env$site == "quarry"], rownames(comm.sum)),])
str(Quarry)
env.q <- q.dat %>% 
  mutate(sev = factor(Burn.Severity, levels=c(1, 2, 3, 4), labels=c("unburned", "low", "mod", "high"), ordered=T), 
         trt = factor(Treatment, levels=c(0, 1), labels=c("control", "treated"))) 

env.q <- env.q[,c("sev", "trt", "NO3", "NH4", "PO4", "plot.cov", "plot.rich", "prop.inv.cov", "prop.ann.cov", "prop.inv.sp", "prop.ann.sp")] %>% 
  set_colnames(c("sev", "trt", "NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover", "prop. annual coverv", "prop. invasive species", "prop. annual species"))

## Stone Canyon 
SC <- as.matrix(comm.sum[match(env$plot[env$site == "sc"], rownames(comm.sum)),])
str(SC)
env.sc <- sc.dat %>% 
  mutate(sev = factor(Burn.Severity, levels=c(1, 2), labels=c("unburned", "low"), ordered=T), 
         trt = factor(Treatment, levels=c(0, 1), labels=c("control", "treated"))) 

env.sc <- env.sc[,c("sev", "trt", "NO3", "NH4", "PO4", "plot.cov", "plot.rich", "prop.inv.cov", "prop.ann.cov", "prop.inv.sp", "prop.ann.sp")] %>% 
  set_colnames(c("sev", "trt", "NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover", "prop. annual coverv", "prop. invasive species", "prop. annual species"))




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


jpeg(here("output/Scree.jpeg"), width = 700, height = 500)
par(mfrow = c(2,2), cex.main = 1.25, cex.lab = 1.25, cex.axis=1.25, cex = 1)
plot(ord.nmds.stress.q,
     ylim = c(0,1),
     ylab = "Stress",
     xlab = "Number of Dimensions",
     main = "Quarry Fire",
     las = 1,
     pch = 16)
abline(h = 0.2)

stressplot(NMDSord.q,
           main = "Quarry")
plot(ord.nmds.stress.sc,
     ylim = c(0,1),
     ylab = "Stress",
     xlab = "Number of Dimensions",
     main = "Stone Canyon",
     las = 1,
     pch = 16)
abline(h = 0.2)

stressplot(NMDSord.sc,
           main = "Stone Canyon")
dev.off()


par(mfrow = c(2,2))




jpeg(here("output/NMDS.jpeg"), width = 700, height = 600)
par(mfrow = c(2,2), cex.main = 1.25, cex.lab = 1.25, cex.axis=1.25, cex=1)

## a - severity
plot(NMDSord.q, type = "n",display = "sites", las = 1) 
points(NMDSord.q, display = "sites", pch = 19, cex = .75, col = pal[as.factor(env.q$sev)])
vectors <- envfit(NMDSord.q, env.q[,c("NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover",  "prop. invasive species")], na.rm = TRUE)
ordiellipse(NMDSord.q, display = "sites", env.q$sev, draw = "lines",
            col = pal, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Unburned", "Low", "Mod", "High"), col = pal, pch = 15, cex = 1, ncol = 2, bty = "n")

## b - fire retardant
plot(NMDSord.q, type = "n",display = "sites", las = 1) 
points(NMDSord.q, display = "sites", pch = 19, cex = .75, col = pal2[as.factor(env.q$trt)])
vectors <- envfit(NMDSord.q, env.q[,c("NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover",  "prop. invasive species")], na.rm = TRUE)
ordiellipse(NMDSord.q, display = "sites", env.q$trt, draw = "lines",
            col = pal2, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Control", "Fire Retardant"), col = pal2, pch = 15, cex = 1, ncol = 2, bty = "n")

## c - severity
plot(NMDSord.sc, type = "n",display = "sites", las = 1) 
points(NMDSord.sc, display = "sites", pch = 19, cex = .75, col = pal[as.factor(env.sc$sev)])
vectors <- envfit(NMDSord.sc, env.sc[,c("NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover",  "prop. invasive species")], na.rm = TRUE)
ordiellipse(NMDSord.sc, display = "sites", env.sc$sev, draw = "lines",
            col = pal, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Unburned", "Low", "Mod", "High"), col = pal, pch = 15, cex = 1, ncol = 2, bty = "n")

## d - fire retardant
plot(NMDSord.sc, type = "n",display = "sites", las = 1) 
points(NMDSord.sc, display = "sites", pch = 19, cex = .75, col = pal2[as.factor(env.sc$trt)])
vectors <- envfit(NMDSord.sc, env.sc[,c("NO3", "NH4", "PO4", "total cover", "total richness", "prop. invasive cover",  "prop. invasive species")], na.rm = TRUE)
ordiellipse(NMDSord.sc, display = "sites", env.sc$trt, draw = "lines",
            col = pal2, label = FALSE) ## ellipses based on slide exposure
plot(vectors, col = "black")
# legend("bottomright", legend = c("Control", "Fire Retardant"), col = pal2, pch = 15, cex = 1, ncol = 2, bty = "n")
dev.off()

