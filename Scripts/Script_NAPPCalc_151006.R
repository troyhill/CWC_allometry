##### This script explores the CWC biomass/NAPP data
# Goals: 
# 1: Plot biomass, stem density, etc. per site per month
# 2: Estimate biomass in mulitple ways
#
# Major change 151006: incorporated un-clipped plots. 
#####



##### load libraries, run previous scripts
#####
# library(ggplot2)
library(reshape2)
library(plyr)
# 'cwc' object has compiled raw allometry data
# run scripts that this code depends on:
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_nondestructiveQuadrats_151008.R")
#####

todaysDate <- substr(as.character(Sys.time()), 1, 10)

 
##### explore data, look for artifacts/errors
#####
hist(cwc$hgt)  # looks reasonable
hist(cwc$mass) 

hist(hgts$hgt) 
hist(hgts$mass)

##### Aggregate data: plot- and marsh-level metrics
#####

### Add quadrat column to CWC and combine with hgts
cwc$quadrat <- "A"

temphgts  <- hgts[, c(1, 2, 4:6, 8:9, 11)]
names(temphgts)[3] <- "monthYear"

# 'quadrats': stem masses from all quadrats
quadrats  <- rbind.fill(list(cwc, temphgts))
tail(quadrats)

# add litter data to standing biomass dataset
lit.proc  <- lit[, c("site", "quadrat", "moYr", "litterMass", "marsh")]
names(lit.proc)[3:4] <- c("monthYear", "mass")
lit.proc$type <- "LITTER"
lit.proc <- seasonLabel(lit.proc, monthYearColumn = "monthYear")
biomass <- rbind.fill(list(quadrats, lit.proc))

CWC.plots <- ddply(biomass, .(marsh, season, monthYear, site, quadrat, type), summarise,
                   mass         =  sum(mass, na.rm = T) / plotSize,
                   stems        =  length(type[!type %in% "LITTER"]) / plotSize,
                   lngth.median =  median(hgt, na.rm = T),
                   lngth.range  =  diff(range(hgt, na.rm = T))
                   ) # warnings come from the infinite values for lngth.range
# head(CWC.plots)
CWC.plots$lngth.range[!is.finite(CWC.plots$lngth.range)] <- NA
CWC.plots$stems[CWC.plots$type %in% "LITTER"] <- NA


### make sure that an absence of LIVE/DEAD/LITTER has not been interpreted as an absence of sampling
### except when there really wasn't any sampling.
### Add absent data and set mass = 0 if one of the three data types is reported (live, dead, litter)
### code prints "missing" observations, so they can be checked. This could impart a severe downward bias to the data, so 
### make sure this is working properly.
nrow(CWC.plots) # 971 151028
for (i in 1:length(unique(CWC.plots$site))) {
  for (j in 1:length(unique(CWC.plots$monthYear))) {
    plot    <- unique(CWC.plots$site)[i]
    time    <- unique(CWC.plots$monthYear)[j]
    subData <- CWC.plots[(CWC.plots$site %in% plot) & (CWC.plots$monthYear == time), ]
    if (nrow(subData) > 0) {
      for (h in 1:length(unique(CWC.plots$quadrat))) {
        subData2 <- subData[subData$quadrat %in% unique(CWC.plots$quadrat)[h], ]
        ### if all three categories are missing for a quadrat, it 
        ### is most likely that the quadrat simply wasn't sampled.
         if (nrow(subData2) == 0) {
           fillData         <- subData[1, ]
           fillData[1, ]    <- NA
           CWC.plots        <- rbind(CWC.plots, fillData[complete.cases(fillData)])
          } else if (nrow(subData2) < 3) {
          # find missing types
          types <- c("LIVE", "DEAD", "LITTER")
          missingTypes <- types[!types %in% subData2$type]
          numberMissing <- length(missingTypes)
          fillData                     <- subData2[1, ]
          if (numberMissing > 1) { # if more than one type is missing, add multiple rows
            for (k in 2:numberMissing) {
              fillData[k, ]            <- subData2[1, ]
            }
          }
          fillData$type                <- missingTypes
          fillData[, grep("mass", names(fillData)):ncol(fillData)]              <- as.numeric(0) # insert NAs from mass column to end of dataset
          fillData[, (grep("mass", names(fillData)) + 2):ncol(fillData)]        <- as.numeric(NA)
          fillData[fillData$type %in% "LITTER", grep("stems", names(fillData))] <- as.numeric(NA)
          CWC.plots <- rbind(CWC.plots, fillData)
          }
        if (exists("fillData")) {
          if (exists("missingData")) {
            missingData <- rbind(missingData, fillData[, c(4, 3, 5, 6, 7)])
          } else if (!exists("missingData")) {
            missingData <- fillData[, c(4, 3, 5, 6, 7)]
          }
          print(fillData[, c("site", "monthYear", "quadrat", "type")])
          rm(fillData)
        }
      }
    }
  }
  rownames(CWC.plots) <- 1:nrow(CWC.plots)
}
nrow(CWC.plots) # 100 larger; 1071
CWC.plots$year <- paste0("20", substr(as.character(CWC.plots$monthYear), 5, 6))
# write.csv(CWC.plots, file = "C:/RDATA/SPAL_allometry/CWC_plotData_151008.csv")
# write.csv(missingData, file = "C:/RDATA/SPAL_allometry/missingData_151008_2.csv") # samples that we are saying has mass = 0
nrow(missingData) # 97
rm(missingData)

### get plot-level metrics (average across quadrats)
cwc.ag <- ddply(CWC.plots[, -c(5)], .(marsh, site, year, season, monthYear, type), colwise(mean, na.rm = T))
cwc.ag.se <- ddply(CWC.plots[, -c(5)], .(marsh, site, year, season, monthYear, type), colwise(se))
head(cwc.ag.se)

cwc.ag <- ddply(CWC.plots, .(marsh, site, year, season, monthYear, type), summarise,
               biomass        =  mean(mass, na.rm = T),
               stemDensity    =  mean(stems, na.rm = T),
               lngth.med      =  mean(lngth.median, na.rm = T),
               lngth.rng      =  mean(lngth.range, na.rm = T),
               biomass.se     =  se(mass),
               stemDensity.se =  se(stems),
               lngth.med.se   =  se(lngth.median),
               lngth.rng.se   =  se(lngth.range)
)

head(cwc.ag)
cwc.ag$time <- as.yearmon(cwc.ag$monthYear, "%b-%y")

# NAs in SE indicate sampling at a single plot? TB-A and TB-B November 2014 (appears correct in the raw data)
# write.csv(cwc.ag, file = "C:/RDATA/SPAL_allometry/CWC_marshData_151008.csv")


### get combined live + dead + litter biomass 
plot.totals <- ddply(CWC.plots, .(marsh, site, quadrat, year, season, monthYear), summarise,
                     biomass.all = sum(mass, na.rm = T),
                     stems.all   = sum(stems, na.rm = T)
)
plot.totals$year <- paste0("20", substr(as.character(plot.totals$monthYear), 5, 6))

# convert from long to wide form
# summary stats across quadrats 
napp <- ddply(CWC.plots, .(marsh, site, year, season, monthYear), summarise,
              live   = mass[type %in% "LIVE"],
              dead   = mass[type %in% "DEAD"],
              litter = mass[type %in% "LITTER"],
              stems.l = stems[type %in% "LIVE"],
              stems.d = stems[type %in% "DEAD"],
              lgth.live = lngth.median[type %in% "LIVE"],
              lgth.dead = lngth.median[type %in% "DEAD"],
              lgth.rng.live = lngth.range[type %in% "LIVE"],
              lgth.rng.dead = lngth.range[type %in% "DEAD"]
              
)
head(napp)

# first, sum standing and fallen dead biomass
totDeadQuads <-  ddply(CWC.plots, .(marsh, site, quadrat, year, monthYear, season), summarise,
                    totalDead = sum(mass[type %in% c("LITTER", "DEAD")]))
# then, average across quadrats
totDeadPlots <-  ddply(totDeadQuads, .(marsh, site, year, monthYear, season), summarise,
                       dead = mean(totalDead, na.rm = T), 
                       dead.se = se(totalDead))
# get other measures on similar basis
nappReplace <- ddply(CWC.plots, .(marsh, site, year, monthYear, season), summarise,
                     live         = mean(mass[type %in% "LIVE"], na.rm = T),
                     deadStems    = mean(mass[type %in% "DEAD"], na.rm = T),
                     litter       = mean(mass[type %in% "LITTER"], na.rm = T),
                     stemDens     =  mean(stems[type %in% "LIVE"], na.rm = T),
                     stemDensDead =  mean(stems[type %in% "DEAD"], na.rm = T),
                lngthMed          =  mean(lngth.median[type %in% "LIVE"], na.rm = T),
                lngthMedDead      =  mean(lngth.median[type %in% "DEAD"], na.rm = T),
                lngthRange        =  mean(lngth.range[type %in% "LIVE"], na.rm = T),
                lngthRangeDead    =  mean(lngth.range[type %in% "DEAD"], na.rm = T),
                live.se           =  se(mass[type %in% "LIVE"]),
                deadStems.se      =  se(mass[type %in% "DEAD"]),
                litter.se         =  se(mass[type %in% "LITTER"]),
#                 dead.se           =  se(mass[type %in% c("DEAD", "LITTER")]) # calculate by sqrt(sum of squared errors)
                stemDens.se       =  se(stems[type %in% "LIVE"]),
                stemDensDead.se   =  se(stems[type %in% "DEAD"]),
                lngthMed.se       =  se(lngth.median[type %in% "LIVE"]),
                lngthMedDead.se   =  se(lngth.median[type %in% "DEAD"]),
                lngthRange.se     =  se(lngth.range[type %in% "LIVE"]),
                lngthRangeDead.se =  se(lngth.range[type %in% "DEAD"])
)
nappReplace$time <- as.yearmon(nappReplace$monthYear, "%b-%y")
tail(nappReplace)
tail(totDeadPlots)

nappReplace$dead    <- totDeadPlots$dead
nappReplace$dead.se <- totDeadPlots$dead.se

napp <- nappReplace

### mean, se by region (for n = 3 sites; could justify using plots as independent replicates)
### combined live + dead + litter biomass, stem density
site.tot <- ddply(plot.totals, .(marsh, monthYear), summarise,
                  biomass    = mean(biomass.all, na.rm = T),
                  stems      = mean(stems.all, na.rm = T),
                  biomass.se = se(biomass.all),
                  stems.se   = se(stems.all)                
)
site.tot$time <- as.yearmon(site.tot$monthYear, "%b-%y")

plot.totals[plot.totals$monthYear %in% "Jun-13", ]
CWC.plots[CWC.plots$monthYear %in% "Jun-13", ]

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
# 
# panel_plot1(y = "biomass", x = "time", ylab = "biomass (g m-2)")
# # ggsave(file = "C:/RDATA/SPAL_allometry/biomass_151008.png", width = wdh, height = hgt, units = "in")
# 
# # september 2015 is extremely high.
# cwc.ag[(cwc.ag$marsh %in% "LUM") & (as.character(cwc.ag$time) %in% "Sep 2015"),]
# # appears to be driven by extremely high live biomass in no-clip plots
# # especially: LUM1-B & -C; LUM3-B & -C
# CWC.plots[(CWC.plots$marsh %in% "LUM") & (as.character(CWC.plots$monthYear) %in% "Sep-15"),]
# hist(biomass$hgt[(biomass$site %in% "LUM1") & (biomass$quadrat %in% "B") & (biomass$monthYear %in% "Sep-15")])
# biomass[(biomass$site %in% "LUM1") & (biomass$quadrat %in% "A") & (biomass$monthYear %in% "Sep-15"), ]
# plot(biomass[(biomass$hgt > 100) & (biomass$type %in% "LIVE"), c("mass", "hgt")])
# plot(cwc[(cwc$hgt > 100) & (cwc$type %in% "LIVE"), c("mass", "hgt")]) # a 103 cm stem weighing 1.9 g?!?
# cwc[(cwc$hgt > 100) & (cwc$type %in% "LIVE"), ]
# 
# 
# 
# 
# panel_plot1(y = "stemDensity", ylab = "stems per square meter")
# # ggsave(file = "C:/RDATA/SPAL_allometry/stemDensity_151008.png", width = wdh, height = hgt, units = "in")
# 
# panel_plot1(y = "lngth.rng", x = "time", ylab = "stem length range")
# # ggsave(file = "C:/RDATA/SPAL_allometry/length3.png", width = wdh, height = hgt, units = "in")
# 
# panel_plot1(y = "lngth.med", x = "time", ylab = "median stem length")
# # ggsave(file = "C:/RDATA/SPAL_allometry/length_median.png", width = wdh, height = hgt, units = "in")


# Only plot LUM, since the other sites are incomplete
m.cwc    <- melt(cwc.ag, id.vars = c("marsh", "site", "time", "year", "season", "monthYear", "type"), 
                 measure.vars = names(cwc.ag)[7:10])
m.cwc.se <- melt(cwc.ag, id.vars = c("marsh", "site", "time", "year", "season", "monthYear", "type"), 
                 measure.vars = names(cwc.ag)[11:14])
m.cwc$value.se <- m.cwc.se$value

lum <- m.cwc[(m.cwc$marsh %in% "LUM"), ]

ggplot(lum, aes(x = as.numeric(time), y = value, col = type)) + geom_point(size = 1.25) + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ site, scale = "free_y") + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())
ggsave(file = paste0("C:/RDATA/SPAL_allometry/LUM_data_", todaysDate, ".png"), width = wdh, height = hgt, units = "in")



y <- "biomass"
x <- "time"
qplot(y = eval(parse(text = y)), x = as.numeric(eval(parse(text = x))), data = site.tot) + geom_point() + 
  geom_errorbar(aes(ymax = eval(parse(text = y)) + eval(parse(text = paste0(y, ".se"))), 
                    ymin = eval(parse(text = y)) - eval(parse(text = paste0(y, ".se")))), width = 0) +
  facet_grid(marsh ~ ., scale='fixed') + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())


m.tot    <- melt(site.tot, id.vars = c("marsh", "time"), measure.vars = names(site.tot)[3:4])
m.tot.se <- melt(site.tot, id.vars = c("marsh", "time"), measure.vars = names(site.tot)[5:6])
m.tot$value.se <- m.tot.se$value

lum.tot <- m.tot[m.tot$marsh %in% "LUM", ]

ggplot(lum.tot, aes(x = as.numeric(time), y = value)) + geom_point() + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ ., scale = "free_y", labeller = labeli) + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())
ggsave(file = paste0("C:/RDATA/SPAL_allometry/LUM_totData_", todaysDate, ".png"), width = wdh, height = hgt, units = "in")

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

# "dead" column in napp is the sum of standing dead and litter
napp2 <- nappCalc(napp, summarize = "TRUE")

napp2$summary
napp2$summary <- marshName(napp2$summary)
plot(napp2$summary[, 5:10])
kable(corstarsl(as.matrix(napp2$summary[, 5:10])))

### this is kind of interesting. Compare with literature values
for (i in 2013:2015) {
  yr <- as.character(i)
  print(kable(corstarsl(as.matrix(napp2$summary[napp2$summary$year == yr, 5:10]))))
  print(noquote(yr))
}


# compare with separate peak standing crop calculation
pscInc <- PSC(napp)
pscInc
napp2$summary[, c(1, 2, 7, 9:10)]

# compare productivity estimates
# first, summarize

dd.napp   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley  = mean(napp.smalley, na.rm = T),
                   MH       = mean(napp.MH, na.rm = T),
                   maxMin   = mean(napp.maxMin, na.rm = T),
                   VTS      = mean(napp.VTS, na.rm = T),
                   psc.live = mean(napp.psc.a, na.rm = T),
                   psc.tot  = mean(napp.psc.b, na.rm = T)
                  )
dd.napp.se   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley.se  = se(napp.smalley),
                   MH.se       = se(napp.MH),
                   maxMin.se   = se(napp.maxMin),
                   VTS.se      = se(napp.VTS),
                   psc.live.se = se(napp.psc.a),
                   psc.tot.se  = se(napp.psc.b)
                   )

m.napp     <- melt(dd.napp, id.vars = c("marsh", "year"))
m.napp.se  <- melt(dd.napp.se, id.vars = c("marsh", "year"))
m.napp$value.se  <- m.napp.se$value


# plots comparing sites across different methods
#####
ggplot(m.napp, aes(x = year, y = value, col = marsh)) + geom_point() + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ ., scale = "fixed", labeller = nappLabelConv) + ylim(0, 4500) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  theme_bw() + theme(legend.title = element_blank())
# ggsave("C:/RDATA/SPAL_allometry/NAPP_compare.png", width = 7, height= 7, units = "in", dpi = 300)

ggplot(m.napp[m.napp$marsh %in% "LUM", ], aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", "Max-Min",
           "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
pngName <- paste0("C:/RDATA/SPAL_allometry/NAPP_compare_LUM_", todaysDate, ".png")
ggsave(pngName, width = 7, height= 4, units = "in", dpi = 300)




ggplot(m.napp, aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ .) +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", "Max-Min",
                                 "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
pngName <- paste0("C:/RDATA/SPAL_allometry/NAPP_compare_allSites_", todaysDate, ".png")
ggsave(pngName, width = 7, height= 4, units = "in", dpi = 300)
#####


#####
##### Explore belowground biomass data (masses are g/m2)
#####

### get max - min for total MOM (similar to Gallagher and Plumley 1979)
### get turnover time (also do this for aboveground biomass): maximum biomass divided by annual increment 
# PSC(bg2, liveCol = "live.bg", deadCol = "dead.bg")
# max total biomass is in 'psc_a' column
bgInc              <- PSC(bg2, liveCol = "total.bg", deadCol = "dead.bg", type = "PSC-A") 
bgInc$turnoverTime <- 12 / (bgInc$psc_a / bgInc$maxMin_a) # max biomass divided by max-min increment. units = months
summary(bgInc)
# MOM turnover time is ~ a month more rapid than aboveground biomass

names(bgInc)[3:5] <- paste0("bg_", names(bgInc)[3:5])
napp2$summary$siteYr <- paste0(pscInc$site, "-", pscInc$year)
bgInc$siteYr <- paste0(bgInc$site, "-", bgInc$year)

abpp <- join_all(list(napp2$summary, pscInc[, c(5:9)], bgInc[, 3:6]), by = "siteYr")
summary(abpp)
abpp <- zeroToNA(abpp)
abpp$rtToSht <- abpp$bg_maxMin_a / abpp$napp.psc.b
summary(abpp[as.character(abpp$year) %in% c("2014"), ])

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

