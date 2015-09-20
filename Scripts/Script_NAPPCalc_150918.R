##### This script explores the CWC data
# Goals: 
# 1: Plot biomass, stem density, etc. per site per month
# 2: Is seasonality at CWC real? 
# 3: Estimate biomass in mulitple ways
# Allometry parameters:
# - Does allometry differ between sites, seasons? 
#     - seasons: plot seasonally-aggregated data together, see if there are differences
#     - sites: use LUM allometry to predict TB biomass, calculate error
# - Compare allometry parameters to literature values
#     - Thursby et al., Morris datasets?
#####



##### load libraries, run previous scripts
#####
# library(ggplot2)
library(reshape2)
library(plyr)
# 'cwc' object has compiled raw allometry data
# source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_mergeData_150918.R")
# or read.csv("C:/RDATA/SPAL_allometry/CWC_allometryData_150918.csv")
#####

##### load ancillary datasets
#####
# litter
lit <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_litter_150917.txt")

# Belowground biomass
bg <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_belowground_150918.txt", skip = 1)
bg <- bg[, c(1, 2, 5, 8, 11, 14)]

# moisture/OM content
om <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_moisture_OM_150917.txt")
om <- om[, c(1:2, 4:6, 8:9, 12)]

# benthic chlorophyll
chl <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_benthicChl_150917.txt", skip = 13)

# nutrients? do these data exist? 

#####

##### declare custom functions
#####
se <- function(x){
  ### Calculates standard errors
  sd(x, na.rm = T) / sqrt(sum(!is.na(x) == T))
}

panel_plot1 <- function(y, x = "time", ylab = "") {
  qplot(y = eval(parse(text = y)), x = as.numeric(eval(parse(text = x))), data = cwc.ag, colour = type) + geom_point() + 
    geom_errorbar(aes(ymax = eval(parse(text = y)) + eval(parse(text = paste0(y, ".se"))), ymin = eval(parse(text = y)) - eval(parse(text = paste0(y, ".se")))), width = 0) +
    facet_grid(marsh ~ ., scale='free_y') + labs(x = "", y = ylab) + 
    theme_bw() + theme(legend.title = element_blank())
}

labeli <- function(variable, value){
  value <- droplevels(value)
  names_li <- list("biomass" = "Biomass (g m-2)", "stemDensity" = "Stem density (m-2)",
                   "length.top3" = "Longest stems (cm)", "stems" = "Stem density (m-2)"
                   )
  return(names_li[value])
}

nappLabelConv <- function(variable, value){
  value <- droplevels(value)
  names_li <- list(
                   "smalley" = "Smalley 1959", "MH" = "Milner Hughes 1968", 
                   "VTS" = "Valiela et al. 1975", 
                   "psc.live" = "Peak (live)",
                   "psc.tot" = "Peak (live + dead)"
  )
  return(names_li[value])
}


#####



##### declare local variables
#####
plotSize <- 0.25^ 2
coreTube <- pi*(6.9/2)^2 # are of coring tube: cm2 
bagMass  <- 5.4 # grams; bags used for belowground cores
  
#####

##### Process litter data
#####
# rename columns
names(lit)    <- c("monthYear", "site", "plot", "quadrat", "plants_tin", "tin", "litterMass")
lit$moYr      <- as.character(lit$monthYear)
lit$monthYear <- as.yearmon(lit$moYr, "%b-%y")
lit$site      <- gsub(" ", "", as.character(lit$site))
lit$plot      <- as.character(lit$plot)
lit$quadrat   <- substr(gsub(" ", "", as.character(lit$quadrat)), 1, 1)

# there are two "A" samples in summer 2014 for LUM1; average them. 
# this removes two rows, but should  only remove one; not sure what the missing row is
lit <- ddply(lit, .(monthYear, site, plot, monthYear, moYr, quadrat), summarise,
             litterMass = mean(litterMass, na.rm = T))
# I'm only interested in quadrat A
lit <- lit[(lit$quadrat %in% "A"), c("monthYear", "moYr", "site", "litterMass")]

#####

##### Process ancillary data
#####
names(bg) <- c("moYr", "site", "depth", "live.bg", "dead.bg", "total.bg")
bg$moYr   <- as.yearmon(bg$moYr, "%b-%y")
bg$site   <- gsub(" ", "", as.character(bg$site))
bg$depth  <- gsub(" ", "", as.character(bg$depth))
bg[, 4:6] <- bg[, 4:6] / coreTube * 10^4 # mass per m2
bg$year <- as.numeric(substr(bg$moYr, 5, 8))
# ignore depth increments for now...
bg2 <- ddply(bg[, c(1:2, 4:6)], .(site, moYr, year), colwise(sum, na.rm = T))
bg2$live.bg[bg2$live.bg == 0] <- NA # zeroes should be NAs (live/dead material not separated)
bg2$dead.bg[bg2$dead.bg == 0] <- NA

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


#####


##### explore data, look for artifacts/errors
#####
hist(cwc$hgt) # verify 150 cm stem? maybe even 130+ cm stems
hist(cwc$mass) # bag scraps are included in this dataset; these are the heavier samples





##### Aggregate data: plot- and marsh-level metrics
#####
CWC.plots <- ddply(cwc, .(site, time, type, marsh, monthYear), summarise,
                   mass         =  sum(mass, na.rm = T) / plotSize,
                   stems        =  length(type) / plotSize,
                   lngth.top3    =  mean(sort(hgt, decreasing = TRUE)[1:3], na.rm = T),
                   lngth.median =  median(hgt, na.rm = T)
                   )

### add litter
head(lit)
head(CWC.plots)

for (i in 1:nrow(lit)) {
  targPlot <- lit$site[i]
  targTime <- lit$moYr[i]
  litterMass <- lit$litterMass[(lit$site %in% targPlot) & (lit$moYr %in% targTime)]
  # add litter mass to dead mass for plot
  CWC.plots$mass[(CWC.plots$type %in% "DEAD") & (CWC.plots$site %in% targPlot) & (CWC.plots$monthYear %in% targTime)] <- litterMass / plotSize + 
    CWC.plots$mass[(CWC.plots$type %in% "DEAD") & (CWC.plots$site %in% targPlot) & (CWC.plots$monthYear %in% targTime)]
}

### make sure that an absence of biomass has not been interpreted as an absence of sampling
nrow(CWC.plots) # 211
for (i in 1:length(unique(CWC.plots$site))) {
  for (j in 1:length(unique(CWC.plots$time))) {
    plot    <- unique(CWC.plots$site)[i]
    time    <- unique(CWC.plots$time)[j]
    subData <- CWC.plots[(CWC.plots$site %in% plot) & (CWC.plots$time == time), ]
    if (nrow(subData) == 1) { 
              fillData                     <- subData[1, ]
              fillData$type                <- ifelse(fillData$type %in% "LIVE", "DEAD", "LIVE")
              fillData[, grep("mass", names(fillData)):ncol(fillData)] <- as.numeric(0) # insert zeroes from mass column to end of dataset
              CWC.plots <- rbind(CWC.plots, fillData)
      }
  }
  rownames(CWC.plots) <- 1:nrow(CWC.plots)
}
nrow(CWC.plots) # 1 larger (212)
CWC.plots$year <- substr(as.character(CWC.plots$time), 5, 9)
# write.csv(CWC.plots, file = "C:/RDATA/SPAL_allometry/CWC_plotData_150918.csv")


### get marsh-level metrics
cwc.ag <- ddply(CWC.plots, .(marsh, time, year, type), summarise,
                   biomass          =  mean(mass, na.rm = T),
                   stemDensity      =  mean(stems, na.rm = T),
                   length.top3      =  mean(lngth.top3, na.rm = T),
                   length.median    =  mean(lngth.median, na.rm = T),
                   biomass.se       =  se(mass),
                   stemDensity.se   =  se(stems),
                   length.top3.se   =  se(lngth.top3),
                   length.median.se =  se(lngth.median)
)

head(cwc.ag)

# NAs in SE indicate sampling at a single plot? TB-A and TB-B November 2014 (appears correct in the raw data)
# write.csv(cwc.ag, file = "C:/RDATA/SPAL_allometry/CWC_marshData_150825.csv")


### get combined live + dead biomass 
plot.totals <- ddply(CWC.plots, .(site, time, year, marsh), summarise,
                     biomass.all = sum(mass, na.rm = T),
                     stems.all   = sum(stems, na.rm = T)
)
plot.totals$year <- substr(as.character(plot.totals$time), 5, 9)



napp <- ddply(CWC.plots, .(marsh, site, time, year), summarise,
              live = mass[type %in% "LIVE"],
              dead = mass[type %in% "DEAD"],
              stems.l = stems[type %in% "LIVE"],
              stems.d = stems[type %in% "DEAD"],
              lgth.live = lngth.top3[type %in% "LIVE"],
              lgth.dead = lngth.top3[type %in% "DEAD"]
)
head(napp)

### mean, se by site
site.tot <- ddply(plot.totals, .(marsh, time), summarise,
                  biomass    = mean(biomass.all, na.rm = T),
                  stems      = mean(stems.all, na.rm = T),
                  biomass.se = se(biomass.all),
                  stems.se   = se(stems.all)                
)



#####


##### Exploratory plots
#####

# our custom function condenses these four lines of code to a single, brief line: 
# qplot(y = biomass, x = as.numeric(time), data = cwc.ag, colour = type) + geom_point() + 
#   geom_errorbar(aes(ymax = biomass + biomass.se, ymin=biomass - biomass.se), width = 0) +
#   facet_grid(marsh ~ ., scale='free_y') + labs(x = "", y = "biomass (g m-2)") + 
#   theme_bw() + theme(legend.title = element_blank())
hgt <- 5
wdh <- 8

panel_plot1(y = "biomass", x = "time", ylab = "biomass (g m-2)")
# ggsave(file = "C:/RDATA/SPAL_allometry/biomass.png", width = wdh, height = hgt, units = "in")


panel_plot1(y = "stemDensity", ylab = "stems per square meter")
# ggsave(file = "C:/RDATA/SPAL_allometry/stemDensity.png", width = wdh, height = hgt, units = "in")

panel_plot1(y = "length.top3", x = "time", ylab = "longest three stems")
# ggsave(file = "C:/RDATA/SPAL_allometry/length3.png", width = wdh, height = hgt, units = "in")

panel_plot1(y = "length.median", x = "time", ylab = "median stem length")
# ggsave(file = "C:/RDATA/SPAL_allometry/length_median.png", width = wdh, height = hgt, units = "in")


# Only plot LUM, since the other sites are incomplete
m.cwc    <- melt(cwc.ag, id.vars = c("marsh", "time", "type"), measure.vars = names(cwc.ag)[5:7])
m.cwc.se <- melt(cwc.ag, id.vars = c("marsh", "time", "type"), measure.vars = names(cwc.ag)[9:11])
m.cwc$value.se <- m.cwc.se[, 5]

lum <- m.cwc[m.cwc$marsh %in% "LUM", ]

ggplot(lum, aes(x = as.numeric(time), y = value, col = type)) + geom_point() + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ ., scale = "free_y", labeller = labeli) + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())
# ggsave(file = "C:/RDATA/SPAL_allometry/LUM_data.png", width = wdh, height = hgt, units = "in")



y <- "biomass"
x <- "time"
qplot(y = eval(parse(text = y)), x = as.numeric(eval(parse(text = x))), data = site.tot) + geom_point() + 
  geom_errorbar(aes(ymax = eval(parse(text = y)) + eval(parse(text = paste0(y, ".se"))), 
                    ymin = eval(parse(text = y)) - eval(parse(text = paste0(y, ".se")))), width = 0) +
  facet_grid(marsh ~ ., scale='free_y') + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())


m.tot    <- melt(site.tot, id.vars = c("marsh", "time"), measure.vars = names(site.tot)[3:4])
m.tot.se <- melt(site.tot, id.vars = c("marsh", "time"), measure.vars = names(site.tot)[5:6])
m.tot$value.se <- m.tot.se$value

lum.tot <- m.tot[m.tot$marsh %in% "LUM", ]

ggplot(lum.tot, aes(x = as.numeric(time), y = value)) + geom_point() + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ ., scale = "free_y", labeller = labeli) + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())
# ggsave(file = "C:/RDATA/SPAL_allometry/LUM_totData.png", width = wdh, height = hgt, units = "in")

#####


##### Compare standing crop estimates
#####
### 1) peak standing crop (max standing biomass (PSC-A: live and PSC-B: live + dead); Linthurst and Reimold 1978 use only live
### 2) Smalley 1959: measure live and dead over time (both +: sum both; 
#                                                   both -: 0; 
#                                                   live + and dead -: live; 
#                                                   live - and dead +: sum both (and if negative, set equal to zero))
# 3) Milner and Hughes 1968 (following Linthurst and Reimold 1978; original text is unclear about the role of dead biomass)
#   sum of positive changes in standing crop over each year
# 4) Valiela, Teal, Sass 1975
#   delB <- change in live standing crop between sampling dates
#   delA <- change in dead standing crop
#   if (delB > 0 and delA < 0) { 
#     e (NAPP increment) <- -delA 
#   } else if (delB < 0) {
#     e <- -(delB + delA)
#    }
#  e shouldn't be negative because assumption is that only live biomass contributes to standing dead material. If that happens, set e to zero
#  if (e < 0) e <- 0


napp2 <- nappCalc(napp, summarize = "TRUE")

napp2$summary
napp2$summary <- marshName(napp2$summary)

# compare with separate peak standing crop calculation
PSC(napp)


# compare productivity estimates
# first, summarize

dd.napp   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley  = mean(napp.smalley, na.rm = T),
                   MH       = mean(napp.MH, na.rm = T),
                   VTS      = mean(napp.VTS, na.rm = T),
                   psc.live = mean(napp.psc.a, na.rm = T),
                   psc.tot  = mean(napp.psc.b, na.rm = T)
                  )
dd.napp.se   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley.se  = se(napp.smalley),
                   MH.se       = se(napp.MH),
                   VTS.se      = se(napp.VTS),
                   psc.live.se = se(napp.psc.a),
                   psc.tot.se  = se(napp.psc.b)
                   )

m.napp     <- melt(dd.napp, id.vars = c("marsh", "year"))
m.napp.se  <- melt(dd.napp.se, id.vars = c("marsh", "year"))
m.napp$value.se  <- m.napp.se$value


# Compare sites across different methods
ggplot(m.napp[!m.napp$year %in% "2015", ], aes(x = year, y = value, col = marsh)) + geom_point() + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ ., scale = "fixed", labeller = nappLabelConv) + ylim(0, 3500) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare.png", width = 7, height= 7, units = "in", dpi = 300)

ggplot(m.napp[(!m.napp$year %in% "2015") & (m.napp$marsh %in% "LUM"), ], aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
           "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare_LUM.png", width = 6, height= 5, units = "in", dpi = 300)




ggplot(m.napp[(!m.napp$year %in% "2015"), ], aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ .) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", 
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare_allSites.png", width = 6, height= 5, units = "in", dpi = 300)



#####
##### Explore belowground biomass data (masses are g/m2)
#####
PSC(bg2, liveCol = "live.bg", deadCol = "dead.bg")


#####


#####
##### Add ancillary data, look for correlations
#####
head(om2)
head(CWC.plots)

 
#####




### TODO:
# run ANOVAs on NAPP methods?
# add ancillary params to napp, napp2$summary

#####
##### 
#####
#####

