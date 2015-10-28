##### This script loads ancillary data & correlates for plant allometry dataset
#####

library(zoo)
library(reshape2)
library(plyr)
library(ggplot2)

##### load ancillary datasets
#####
# litter
lit <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_litter_151020.txt")

# Belowground biomass
bg <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_belowground_150922.txt", skip = 1)
bg <- bg[, c(1, 2, 5, 8, 11, 14)]

# moisture/OM content
om <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_moisture_OM_150917.txt")
om <- om[, c(1:2, 4:6, 8:9, 12)]

# benthic chlorophyll
chl <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_benthicChl_150925.txt", skip = 13,
                  na.strings = "#DIV/0!")
chl <- chl[!is.na(chl[, 18]), c(1:2, 5, 18:19, 22:23)] # leaving out irrelevant data (meter mark, quadrat)

# nutrients 
nutrientDirectory <- "C:/RDATA/SPAL_allometry/data_LUM123/nutrients/"
fileList <- list.files(nutrientDirectory)
nutFiles <- paste0(nutrientDirectory, fileList[grep("nutrients", list.files(nutrientDirectory))])
test <- lapply(nutFiles, read.delim, skip = 1, na.strings = "#DIV/0!")
head(test[[5]])
nuts <- do.call("rbind", test)



#####

chlSyringeVol <- chl_syringeArea * 0.5


##### Process litter data
#####
# rename columns
names(lit)    <- c("monthYear", "site", "plot", "quadrat", "plants_tin", "tin", "litterMass")
lit$moYr      <- as.character(lit$monthYear)
lit$monthYear <- as.yearmon(lit$moYr, "%b-%y")
lit$site      <- gsub(" ", "", as.character(lit$site))
lit$plot      <- as.character(lit$plot)
lit$quadrat   <- substr(gsub(" ", "", as.character(lit$quadrat)), 1, 1)


### add bag scraps (in cwc) to quadrat A's litter (lit)
head(bagScraps)
names(bagScraps)  <- c("site", "monthYear", "ID", "type", "hgt", "litterMass", "stemMass", "leafMass", "marsh", "moYr")
bagScraps$quadrat <- "A"
bagScraps$site <- as.character(bagScraps$site)

lit <- rbind(lit[, which(names(lit) %in% names(bagScraps))], bagScraps[, which(names(bagScraps) %in% names(lit))])

# average duplicated samples within a quadrat (e.g., from bag scraps and litter) 
lit <- ddply(lit, .(site, quadrat, monthYear, moYr), summarise,
             litterMass = mean(litterMass, na.rm = T))
lit <- marshName(lit)

#####

##### Process ancillary data
#####

### belowground biomass data
names(bg) <- c("moYr", "site", "depth", "live.bg", "dead.bg", "total.bg")
bg$moYr   <- as.yearmon(bg$moYr, "%b-%y")
bg$site   <- gsub(" ", "", as.character(bg$site))
bg$depth  <- gsub(" ", "", as.character(bg$depth))
bg[, 4:6] <- bg[, 4:6] / coreTube * 10^4 # mass per m2
bg$year   <- as.numeric(substr(as.character(bg$moYr), 5, 8))
# remove suspicious LUM2 samples from Jan 2014
bg <- bg[-c(grep("LUM2", bg$site)[grep("LUM2", bg$site) %in% grep("Jan 2014", as.character(bg$moYr))]), ]

# ignore depth increments for now...
bg2 <- ddply(bg[, c(1:2, 4:7)], .(site, moYr, year), colwise(sum, na.rm = T))
bg2$year   <- as.numeric(substr(as.character(bg2$moYr), 5, 8))
bg2$live.bg[bg2$live.bg == 0] <- NA # zeroes should be NAs (live/dead material not separated)
bg2$dead.bg[bg2$dead.bg == 0] <- NA


### organic matter, water content data
names(om) <- c("moYr", "site", "depth", "bagSoilWet", "tin", "tinPlusWetSample", "SoilTin70C", "ashedSoilTin")
om$moYr <- as.yearmon(om$moYr, "%b-%y")
om$site <- gsub(" ", "", as.character(om$site))
# note: moisture content is calculated wrong in spreadsheet
om$gravWtrCont  <- (om$tinPlusWetSample - om$SoilTin70C) / (om$tinPlusWetSample - om$tin)
om$volWtrCont   <- (om$tinPlusWetSample - om$SoilTin70C) / (coreTube * 5) # assumes 1.0 g/cm3 water density
om$pctOM        <- (om$SoilTin70C - om$ashedSoilTin) / (om$SoilTin70C - om$tin)
om$dryMass      <- (om$bagSoilWet - bagMass) * (1-om$gravWtrCont) # dry mass (g; 70C)
om$bulkDens     <- om$dryMass / (coreTube*5)
om$OMVol        <- (om$pctOM * om$dryMass)/1.1
om$IMVol        <- ((1-om$pctOM) * om$dryMass)/2.65
om$wtrVol       <- om$volWtrCont*(coreTube * 5)
om$poreVol      <- (coreTube * 5) - om$wtrVol - om$IMVol - om$OMVol
om$pctPoreSpace <- om$poreVol / (coreTube * 5)
om2 <- om[, c("site", "moYr", "pctOM", "bulkDens", "volWtrCont", "poreVol", "wtrVol", "IMVol", "OMVol")]


### chlorophyll A & phaeopigment inventories, concentrations
names(chl)       <- c("date", "plot", "sampleName", "chlA_inv", "pgmt_inv", "chlA_conc", "pgmt_conc") # chlorophyll A and pigments, expressed in per cm2 and per gram dry mass of sediment
chl$chlA_vol     <- (chl$chlA_inv * chl_syringeArea) / chlSyringeVol # grams per cm3 
chl$pgmt_vol     <- (chl$chlA_inv * chl_syringeArea) / chlSyringeVol
chl$chlA_conc    <- chl$chlA_conc / 1e3 # convert concentrations from ug/g to mg/g
chl$pgmt_conc    <- chl$pgmt_conc / 1e3
chl$chlA_inv     <- chl$chlA_inv / 1e3 * 1e4 # convert inventories from ug/cm2 to mg/cm2
chl$pgmt_inv     <- chl$pgmt_inv / 1e3 * 1e4
chl$plot         <- gsub(" ", "", as.character(chl$plot))
chl$sampleName   <- gsub(" ", "", as.character(chl$sampleName))
chl$moYr         <- as.yearmon(chl$date, "%y-%b")
chl              <- marshName(chl, siteCol = "plot")
# average replicate samples, where available
chl2    <- ddply(chl[, c(2, 4:ncol(chl))], .(plot, moYr, marsh), colwise(mean, na.rm = T))
chl2.se <- ddply(chl[, c(2, 4:ncol(chl))], .(plot, moYr, marsh), colwise(se))
names(chl2.se) <- paste0(names(chl2.se), ".se")
chl2    <- cbind(chl2, chl2.se[, 4:7])

### nutrients
names(nuts) <- c("ID", "po4", "no3", "si", "nh4")
# extract plot name (called "site" for consistency)
int         <- strsplit(as.character(nuts$ID), "-")
nuts$site   <- gsub(" ", "", sapply(int, "[", 1)) # extracts first element from each vector in the list
nuts        <- marshName(nuts)

# extract date
indices   <- regexpr("/", as.character(nuts$ID))
int       <- gsub(" ", "", substr(as.character(nuts$ID), indices - 2, indices + 2))
mos       <- month.abb[as.numeric(sapply(strsplit(int, "/"), "[", 1))]
yr        <- paste0("20", sapply(strsplit(int, "/"), "[", 2))
nuts$moYr <- as.yearmon(paste0(mos, "-", yr), "%b-%Y")


# how to do this with sapply and strsplit?
# labelVec   <- c("site1 BAY DECOMP PROJ 11/14", "site1 DECOMP 6/14", "site1 BAY DECOMP 5/14")
# labelList  <- strsplit(labelVec, " ")
# lengths      <- sapply(labelList, FUN = length)
# sapply(labelList, "[", c(lengths)) 


# extract location, if stated
locIndices   <- regexpr("-", as.character(nuts$ID))
nuts$loc     <- gsub(" ", "", substr(as.character(nuts$ID), locIndices + 1, locIndices + 3))
nuts$loc     <- ifelse(nuts$loc %in% c("5", "10", "20"), paste0(nuts$loc, "M"), nuts$loc)
# change "10M" to "20M" for sake of comparisons
nuts$loc     <- ifelse(nuts$loc %in% c("10M"), "20M", nuts$loc)


# write.csv(nuts, file = "C:/RDATA/plantsNutrientData.csv")
#####




#####
##### Plots exploring ancillary data
#####

### nutrient data
### analysis including all plots (all "marshes")
m.nuts <- melt(nuts[nuts$loc %in% c("1B", "4B"), 2:ncol(nuts)], id.vars = c("site", "moYr", "loc", "marsh"))

ggplot(m.nuts, aes(x = as.numeric(moYr), y = value, col = site)) + geom_point() + 
  facet_grid(variable ~ loc, scales = "free_y") + theme_bw() + 
  labs(y = expression("Concentration ("*mu*"M)"), x = "")

### focus on Bay samples. Aggregate to region-scale and plot with error bars
nut.prep <- ddply(nuts[nuts$loc %in% c("BAY"), c(2:5, 7, 8)], .(marsh, moYr), colwise(mean, na.rm = T))
nut.se   <- ddply(nuts[nuts$loc %in% c("BAY"), c(2:5, 7, 8)], .(marsh, moYr), colwise(se))
m.bay    <- melt(nut.prep, id.vars = c("moYr", "marsh"))
m.bay.se <- melt(nut.se,   id.vars = c("moYr", "marsh"))
m.bay$se <- m.bay.se$value

ggplot(m.bay, aes(x = as.numeric(moYr), y = value, col = marsh)) + geom_point() + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0) +
  facet_grid(variable ~ ., scales = "free_y") + theme_bw() + 
  labs(y = expression("Concentration ("*mu*"M)"), x = "")
# ggsave("C:/RDATA/SPAL_allometry/BayNutrients_150925.png", width = 6, height= 5, units = "in", dpi = 300)



### benthic chlorophyll, by plot
m.chl    <- melt(chl2[, 1:7], id.vars = c(1:3))
m.chl.se <- melt(chl2[, c(1:3, 8:11)], id.vars = c(1:3))
m.chl$se <- m.chl.se$value

ggplot(m.chl, aes(x = as.numeric(moYr), y = value, col = plot)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0) +
  facet_grid(variable ~ ., scales = "free_y", labeller = chlLabeller) + theme_bw() + 
  labs(y = "Benthic chlorophyll", x = "")
# ggsave("C:/RDATA/SPAL_allometry/chlA_150925.png", width = 5, height= 7, units = "in", dpi = 300)







### OM, water content, bulk density


#####
