##### This script explores the CWC biomass/NAPP data
# Goals: 
# 1: Plot biomass, stem density, etc. per site per month
# 2: Estimate biomass in mulitple ways
#
#####



##### load libraries, run previous scripts
#####
library(corrgram)
library(corrplot)
library(reshape2)
library(plyr)
library(dplyr) # for bind_rows()
library(Hmisc) # for function rcorr()
# 'cwc' object has compiled raw allometry data
# run scripts that this code depends on:
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_nondestructiveQuadrats_151008.R")
par(fig = c(0,1,0,1), mar = c(4, 4.5, 0.5, 0.5))
#####

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

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

# columns to remove
hgt_cols <- which(!names(hgts) %in% c("seas", "notes", "analyst", "mass", "mass2", "mass3"))

### which mass column to use (differ by allometry equations used to estimate masses)
mass_col <- which(names(hgts) %in% c("mass2"))

temphgts  <- hgts[, c(hgt_cols, mass_col)]
names(temphgts)[3] <- "monthYear"
names(temphgts)[ncol(temphgts)] <- "mass"

# 'quadrats': stem masses from all quadrats
quadrats  <- rbind.fill(list(cwc, temphgts))


# add litter data to standing biomass dataset
lit.proc  <- lit[, c("site", "quadrat", "moYr", "litterMass", "marsh")]
names(lit.proc)[3:4] <- c("monthYear", "mass")
lit.proc$type <- "LITTER"
lit.proc <- seasonLabel(lit.proc, monthYearColumn = "monthYear")
biomass <- rbind.fill(list(quadrats, lit.proc))
biomass <- seasonLabel(biomass, monthYearColumn = "monthYear")
biomass$seas  <- substr(biomass$season, 1, 4)

suppressWarnings(CWC.plots <- ddply(biomass, .(marsh, season, monthYear, site, quadrat, type), summarise,
                   mass         =  sum(mass, na.rm = T) / plotSize,
                   stems        =  length(type[!type %in% "LITTER"]) / plotSize,
                   lngth.median =  median(hgt, na.rm = T),
                   lngth.range  =  diff(range(hgt, na.rm = T))
                   ) # warnings come from the infinite values for lngth.range
)

CWC.plots$lngth.range[!is.finite(CWC.plots$lngth.range)] <- NA
CWC.plots$stems[CWC.plots$type %in% "LITTER"] <- NA
# odd points from CWC.plots (there are many...) 
# 789  TB-A sumr 13    Jun-13  TB2       B   LIVE    0.000   368        68.00        53.0 2013 sumr
# 790  TB-A sumr 13    Jun-13  TB2       C   DEAD    0.000   128        30.50        52.5 2013 sumr
# 791  TB-A sumr 13    Jun-13  TB2       C LITTER   30.400    NA           NA          NA 2013 sumr
# biomass[(biomass$site %in% "TB2") & (biomass$type %in% "LIVE") & (biomass$monthYear %in% "Jun-13"), ]

### make sure that an absence of LIVE/DEAD/LITTER has not been interpreted as an absence of sampling
### except when there really wasn't any sampling.
### Add absent data and set mass = 0 if one of the three data types is reported (live, dead, litter)
### code prints "missing" observations, so they can be checked. This could impart a severe downward bias to the data, so 
### make sure this is working properly.
nrow(CWC.plots) #  151103
for (i in 1:length(unique(CWC.plots$site))) {
  for (j in 1:length(unique(CWC.plots$monthYear))) {
    plot    <- unique(CWC.plots$site)[i]
    time    <- unique(CWC.plots$monthYear)[j]
    subData <- CWC.plots[(CWC.plots$site %in% plot) & (CWC.plots$monthYear == time), ]
    if (nrow(subData) > 0) {
      for (h in 1:length(unique(CWC.plots$quadrat))) {
        subData2 <- subData[subData$quadrat %in% unique(CWC.plots$quadrat)[h], ]
        ### Assumption: if all three categories are missing for a quadrat, it 
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
nrow(CWC.plots) 
CWC.plots$year <- paste0("20", substr(as.character(CWC.plots$monthYear), 5, 6))
# write.csv(CWC.plots, file = "C:/RDATA/SPAL_allometry/CWC_plotData_151127.csv")
# write.csv(missingData, file = "C:/RDATA/SPAL_allometry/missingData_151008_2.csv") # samples that we are saying has mass = 0
nrow(missingData) # 41: 160111
### May 2014 - no litter, anywhere?
rm(missingData)


CWC.plots       <- seasonLabel(CWC.plots, monthYearColumn = "monthYear")
CWC.plots$seas  <- substr(CWC.plots$season, 1, 4)

### get plot-level metrics (average across quadrats)
colsToRemove <- which(names(CWC.plots) %in% c("quadrat", "seas"))
cwc.ag <- ddply(CWC.plots[, -c(colsToRemove)], .(marsh, site, year, season, monthYear, type), colwise(mean, na.rm = T))
cwc.ag.se <- ddply(CWC.plots[, -c(colsToRemove)], .(marsh, site, year, season, monthYear, type), colwise(se))

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

# head(cwc.ag)
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
# test <- createUniqueID(CWC.plots, c("site", "type", "monthYear"))
# head(test)
# a <- unique(test[duplicated(test)])
# test[test %in% a[1]]

napp <- ddply(CWC.plots, .(marsh, site, year, season, monthYear), summarise,
              live   = mean(mass[type %in% "LIVE"], na.rm = TRUE),
              dead   = mean(mass[type %in% "DEAD"], na.rm = TRUE),
              litter = mean(mass[type %in% "LITTER"], na.rm = TRUE),
              stems.l = mean(stems[type %in% "LIVE"], na.rm = TRUE),
              stems.d = mean(stems[type %in% "DEAD"], na.rm = TRUE),
              lgth.live = mean(lngth.median[type %in% "LIVE"], na.rm = TRUE),
              lgth.dead = mean(lngth.median[type %in% "DEAD"], na.rm = TRUE),
              lgth.rng.live = mean(lngth.range[type %in% "LIVE"], na.rm = TRUE),
              lgth.rng.dead = mean(lngth.range[type %in% "DEAD"], na.rm = TRUE)
              
)
# head(napp)

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
# tail(nappReplace)
# tail(totDeadPlots)

nappReplace$dead    <- totDeadPlots$dead
nappReplace$dead.se <- totDeadPlots$dead.se

napp <- nappReplace
napp <- napp[order(napp$marsh, napp$site, napp$time), ]

### mean, se by region (for n = 3 sites; could justify using plots as independent replicates)
### combined live + dead + litter biomass, stem density
site.tot <- ddply(plot.totals, .(marsh, monthYear), summarise,
                  biomass    = mean(biomass.all, na.rm = T),
                  stems      = mean(stems.all, na.rm = T),
                  biomass.se = se(biomass.all),
                  stems.se   = se(stems.all)                
)
site.tot$time <- as.yearmon(site.tot$monthYear, "%b-%y")

# plot.totals[plot.totals$monthYear %in% "Jun-13", ]
# CWC.plots[CWC.plots$monthYear %in% "Jun-13", ]

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



##### Combine and explore above and belowground datasets
#####
cwc.ag2 <- cwc.ag

bg2$ID <- createUniqueID(bg2, c("site", "moYr"))
napp$ID <- createUniqueID(napp, c("site", "time"))
cwc.ag2$ID <- createUniqueID(cwc.ag2, c("site", "time"))
cwc.ag2$year <- as.numeric(cwc.ag2$year)
combd <- join_all(list(napp, bg2[, -c(1:3)]), by = "ID") # should maybe combine with cwc.ag or CWC.plots?? 

bg2_tmp <- bg2[, c(1:3, 6:7)]
bg2_tmp$type <- "BELOWGRD"
bg2_tmp$time <- as.numeric(bg2_tmp$moYr)
names(bg2_tmp)[4] <- "biomass"
comb <- bind_rows(list(cwc.ag2, bg2_tmp))
comb <- marshName(comb[, -1])
comb$numTime <- as.numeric(comb$time)

comb.ag <- ddply(comb[comb$type %in% "BELOWGRD", ], .(marsh, year, numTime), summarise,
                 MOM        =  mean(biomass, na.rm = T),
                 MOM.se     =  se(biomass))


m.comb    <- melt(comb, id.vars = c("marsh", "site", "year", "time", "type"), 
                 measure.vars = names(comb)[6:9])
m.comb.se <- melt(comb, id.vars = c("marsh", "site", "year", "time", "type"), 
                 measure.vars = names(comb)[10:13])
m.comb$value.se <- m.comb.se$value

comb.lum <- m.comb[(m.comb$marsh %in% "LUM"), ]
comb.lum$type <- ordered(comb.lum$type, c("LIVE", "DEAD", "LITTER", "BELOWGRD"))



comb.lum.tot <- comb.lum[comb.lum$type %in% "BELOWGRD", ]
comb.lum.tot$variable <- "belowground" 
comb.lum.tot$type <- "total"
comb.lum.tot$value.se <- comb.lum.tot$value * 0.14 # 14% variation from replicated cores in 2013
# make live dataset
comb.lum.live <- comb.lum.tot
comb.lum.live$type <- "LIVE"
comb.lum.live$value <- comb.lum.live$value * 0.80
comb.lum.live$value.se <- comb.lum.live$value * 0.14

# make dead dataset
comb.lum.dead <- comb.lum.tot
comb.lum.dead$type <- "DEAD"
comb.lum.dead$value <- comb.lum.dead$value * 0.20
comb.lum.dead$value.se <- comb.lum.dead$value * 0.14


# merge datasets and plot
lumTemp <- lum[!lum$variable %in% "lngth.rng", ]
lumTemp2 <- rbind.fill(list(lumTemp, comb.lum.live, comb.lum.dead))

lumTemp2$variable <- ordered(lumTemp2$variable, c("biomass", "stemDensity", "lngth.med", "belowground"))

lumTemp2$variable <- revalue(lumTemp2$variable, c("biomass" = "Aboveground ~(g %.% m ^2)", 
                             "stemDensity" = "Stems ~(stems %.% m ^2)",
                             "lngth.med"   = "Median~length~(cm)",
                             "belowground"   = "Belowground~(g %.% m ^-2)"))

lumTemp2$site <- revalue(lumTemp2$site, c("LUM1" = "LUM1",
"LUM2" = "LUM2",
"LUM3" = "LUM3"))

ggplot(lumTemp2[!lumTemp2$variable %in% "Median~length~(cm)", ], aes(x = as.numeric(time), y = value, col = type)) + geom_point(size = 1.25) + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se), width = 0) +
  facet_grid(variable ~ site, scale = "free_y", labeller = label_parsed) + labs(x = "", y = "") + 
  theme_bw() + theme(legend.title = element_blank())

ggsave(file = paste0("C:/RDATA/SPAL_allometry/LUM_biomassTrends_3rows", todaysDate, ".png"), width = 7, height = 7, units = "in")







bg2 <- marshName(bg2)
y <- bg2$live.bg[!is.na(bg2$dead.bg)] / 1e3
x <- bg2$total.bg[!is.na(bg2$dead.bg)] / 1e3
summary(lm1 <- lm(y ~ x))
plot(y ~ x, 
     ylab = expression("live biomass (kg"%.%"m"^-2~")"),
     xlab = expression("total biomass (kg"%.%"m"^-2~")"),
     pch = 19, cex = 0.7, las = 1
)
points(y = bg2$live.bg[bg2$marsh %in% "LUM"] / 1e3, x = bg2$total.bg[bg2$marsh %in% "LUM"] / 1e3, pch = 19, col = "red")
points(y = bg2$live.bg[bg2$marsh %in% "TB-A"] / 1e3, x = bg2$total.bg[bg2$marsh %in% "TB-A"] / 1e3, pch = 19, col = "green")
abline(lm1, lty = 2)
text(x = 3, y = 11, pos = 4, "LUM", col = "red")
text(x = 3, y = 10, pos = 4, "TB-A", col = "green")
text(x = 3, y = 9, pos = 4, "TB-B")

summary((bg2$live.bg / 1e3) / (bg2$total.bg / 1e3))
summary((bg2$live.bg / 1e3) / (bg2$dead.bg / 1e3))


plot(I(MOM / 1e3) ~ numTime, data = comb.ag[comb.ag$marsh %in% "LUM",],
     ylab = expression("total belowground biomass (kg"%.%"m"^-2~")"),
     xlab = expression(""),
     pch = 19, cex = 0.7, las = 1,
     ylim = c(0, 10)
)
segments(comb.ag$numTime[comb.ag$marsh %in% "LUM"], I((comb.ag$MOM[comb.ag$marsh %in% "LUM"] + comb.ag$MOM.se[comb.ag$marsh %in% "LUM"]) / 1e3), comb.ag$numTime[comb.ag$marsh %in% "LUM"], I((comb.ag$MOM[comb.ag$marsh %in% "LUM"] - comb.ag$MOM.se[comb.ag$marsh %in% "LUM"]) / 1e3))




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
napp2 <- nappCalc(napp, summarize = "TRUE", EOSL = "TRUE", EOSL_window = 0)

# napp2$summary
napp2$summary <- marshName(napp2$summary)
plot(napp2$summary[, c(5:10, 17)])
plot(napp2$summary[napp2$summary$marsh %in% "LUM", c(5:10, 17)])
# correlations for LUM (all years)
kable(corstarsl(as.matrix(napp2$summary[napp2$summary$marsh %in% "LUM", c(5:6, 8:9, 17)]))) # n = 8-21 (each site-year counted as independent obs)

corrgram(napp2$summary[, c(5:10, 17)], order = NULL, lower.panel = panel.shade, 
         upper.panel = NULL, text.panel = panel.txt, 
         main = "")



M    <- cor(napp2$summary[, c(5:6, 8:9, 17)], use = "complete.obs")
res1 <- cor.mtest(napp2$summary[, c(5:6, 8:9, 17)], 0.95)
corrplot(M, p.mat = res1[[1]], method = "circle", type = "lower", insig = "blank")

kable(corstarsl(as.matrix(napp2$summary[, c(5:6, 8:9, 17)])))



# correlations for all sites (2013 & 2014)
kable(corstarsl(as.matrix(napp2$summary[((napp2$summary$year != 2015)), c(5:10, 17)])))



### this is kind of interesting (except n = 2 for some comparisons). Compare with literature values
for (i in 2013:2015) {
  yr <- as.character(i)
  print(kable(corstarsl(as.matrix(napp2$summary[(napp2$summary$year == yr), c(5:10, 17)]))))
  print(noquote(yr))
}


# compare with separate peak standing crop calculation
pscInc <- PSC(napp)
pscInc[1:4, ]
napp2$summary[1:4, c(1, 2, 7, 9:10)]

# compare productivity estimates
# first, summarize

dd.napp   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley  = mean(napp.smalley, na.rm = T),
                   MH       = mean(napp.MH, na.rm = T),
                   maxMin   = mean(napp.maxMin, na.rm = T),
                   VTS      = mean(napp.VTS, na.rm = T),
                   psc.live = mean(napp.psc.a, na.rm = T),
                   psc.tot  = mean(napp.psc.b, na.rm = T),
                   eosl     = mean(napp.EOSL, na.rm = T)
                  )
dd.napp.se   <- ddply(napp2$summary, .(marsh, year), summarise,
                   smalley.se  = se(napp.smalley),
                   MH.se       = se(napp.MH),
                   maxMin.se   = se(napp.maxMin),
                   VTS.se      = se(napp.VTS),
                   psc.live.se = se(napp.psc.a),
                   psc.tot.se  = se(napp.psc.b),
                   eosl.se     = se(napp.EOSL)
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

ggplot(m.napp[(m.napp$marsh %in% "LUM") & (m.napp$variable %in% unique(m.napp$variable)[c(1:2, 4:5, 7)]), ], aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", "Valiela et al. 1975", 
                                 "Peak (live)", "End-of-season live")) +
  theme_bw() + theme(legend.title = element_blank())
pngName <- paste0("C:/RDATA/SPAL_allometry/NAPP_compare_LUM_", todaysDate, ".png")
ggsave(pngName, width = 7, height= 4, units = "in", dpi = 300)




ggplot(m.napp[m.napp$variable %in% unique(m.napp$variable)[c(1:2, 4:5, 7)], ], 
       aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(marsh ~ ., scale = "free_y") +
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
  scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", "Valiela et al. 1975", 
                                 "Peak (live)", "End-of-season live")) +
  theme_bw() + theme(legend.title = element_blank())
pngName <- paste0("C:/RDATA/SPAL_allometry/NAPP_compare_allSites_", todaysDate, ".png")
ggsave(pngName, width = 6, height= 6, units = "in", dpi = 300)




ggplot(m.napp[(m.napp$variable %in% c("smalley", "VTS", "eosl")) & (m.napp$marsh %in% "LUM"),], aes(x = year, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = value - value.se, ymax = value + value.se),
                width = 0, position = position_dodge(width = 0.9)) + 
  labs(x = "", y = expression("NAPP (g "%.%m^-2%.%yr^-1~")")) + 
    scale_fill_discrete(labels = c("Smalley 1959", 
                                   "Valiela et al. 1975", "End-of-season live")) +
  #   scale_fill_discrete(labels = c("Smalley 1959", "Milner & Hughes 1968", "Max-Min",
#                                  "Valiela et al. 1975", "Peak (live)", "Peak (live + dead)")) +
  theme_bw() + theme(legend.title = element_blank())
ggsave("C:/RDATA/SPAL_allometry/NAPP_3_LUM_only.png", width = 7, height= 4, units = "in", dpi = 300)

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
# summary(bgInc)
# MOM turnover time is ~ a month more rapid than aboveground biomass

names(bgInc)[3:5] <- paste0("bg_", names(bgInc)[3:5])
napp2$summary$siteYr <- paste0(pscInc$site, "-", pscInc$year)
bgInc$siteYr <- paste0(bgInc$site, "-", bgInc$year)

abpp <- join_all(list(napp2$summary, bgInc[, 3:6]), by = "siteYr")
abpp <- zeroToNA(abpp)
abpp$year <- as.numeric(abpp$year)
abpp$rtToSht <- abpp$bg_maxMin_a / abpp$napp.psc.b # ratio of below:aboveground productivity (0.5 - 6x)
# summary(abpp)


qplot(y = rtToSht, x = year, data = abpp, col = site)
ddply(abpp[, c(2:11, 17, 18, 20:21, 23)], .(marsh, year), numcolwise(mean))
ddply(abpp[, c(2:11, 17, 18, 20:21, 23)], .(marsh, year), numcolwise(se))


#####


#####
##### Add ancillary data, look for correlations
#####

 
#####




### TODO:
# run ANOVAs on NAPP methods?
# add ancillary params to napp, napp2$summary




### Create GRIID-C output 
### individual stem masses & heights
quadNames <- which(names(quadrats) %in% c("marsh", "site", "time", "quadrat", "type", "hgt", "mass", "ID"))
gridC_allom <- quadrats[, quadNames]
gridC_allom$year  <- substr(gridC_allom$time, 5, 9)
gridC_allom$month <- substr(gridC_allom$time, 1, 3)
# latitude
gridC_allom$latitude[gridC_allom$site %in% "TB1"]  <- 29.29316111
gridC_allom$latitude[gridC_allom$site %in% "TB2"]  <- 29.29139722
gridC_allom$latitude[gridC_allom$site %in% "TB3"]  <- 29.30445
gridC_allom$latitude[gridC_allom$site %in% "TB4"]  <- 29.29659167
gridC_allom$latitude[gridC_allom$site %in% "LUM1"] <- 29.2554
gridC_allom$latitude[gridC_allom$site %in% "LUM2"] <- 29.25605278
gridC_allom$latitude[gridC_allom$site %in% "LUM3"] <- 29.25796944
  # longitude
gridC_allom$longitude[gridC_allom$site %in% "TB1"]  <- -90.60505
gridC_allom$longitude[gridC_allom$site %in% "TB2"]  <- -90.60300556
gridC_allom$longitude[gridC_allom$site %in% "TB3"]  <- -90.56475833
gridC_allom$longitude[gridC_allom$site %in% "TB4"]  <- -90.54983333
gridC_allom$longitude[gridC_allom$site %in% "LUM1"] <- -90.66478889
gridC_allom$longitude[gridC_allom$site %in% "LUM2"] <- -96.818075
gridC_allom$longitude[gridC_allom$site %in% "LUM3"] <- -90.66144444

gridC_allom$dist  <- 20
gridC_allom$dist[gridC_allom$site %in% "LUM3"] <- 10


gridC_allom$marsh[gridC_allom$marsh %in% "TB-B"] <- "Lake Barre"
gridC_allom$marsh[gridC_allom$marsh %in% "TB-A"] <- "Bay LaFleur"
gridC_allom$marsh[gridC_allom$marsh %in% "LUM"] <- "LUMCON"
gridC_allom$mass[gridC_allom$quadrat %in% c("B", "C")] <- NA
gridC_allom$ID <- toupper(gridC_allom$ID)

gridC_allom <- gridC_allom[order(gridC_allom$marsh, gridC_allom$site, gridC_allom$time, gridC_allom$type, gridC_allom$quadrat), ]

dataset11 <- gridC_allom[, c("year", "month", "marsh", "site", "latitude", "longitude", "dist", "quadrat", "ID", "type", "hgt", "mass")]

# write.csv(dataset11[dataset11$year < 2015, ], file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC1/dataset11_160105.csv")
# write.csv(dataset11[dataset11$year == 2015, ], file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC2/dataset11_160105.csv")


### allometry parameters
seasonRegionParams
# write.csv(seasonRegionParams, file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC1/AllometryParameters_151213.csv")

### aggregate quadrat-level data
uniqueIDcols <- c("site", "quadrat", "monthYear", "dist")
quadData <- ddply(CWC.plots, .(year, monthYear, marsh, site, quadrat), summarise,
                  # first, live-stem data
                  stemDens.l       =  stems[type %in% "LIVE"],
                  lngth.med.l      =  lngth.median[type %in% "LIVE"],
                  biomass.l        =  mass[type %in% "LIVE"], 
                  stemDens.d       =  stems[type %in% "DEAD"],
                  lngth.med.d      =  lngth.median[type %in% "DEAD"],
                  biomass.d        =  mass[type %in% "DEAD"],
                  biomass.litter   =  mass[type %in% "LITTER"]
)
# add latitude & longitude
for (i in 1:nrow(quadData)) {
  quadData$latitude[i]  <- gridC_allom$latitude[gridC_allom$site %in% quadData$site[i]][1]
  quadData$longitude[i] <- gridC_allom$longitude[gridC_allom$site %in% quadData$site[i]][1]
}

quadData$monthYear <- as.character(as.yearmon(quadData$monthYear, "%b-%y"))
quadData$dist <- 20
quadData$dist[quadData$site %in% "LUM3"] <- 10
quadData$ID   <- createUniqueID(quadData, uniqueIDcols)

# belowground biomass
bg.quad  <- zeroToNA(bg_int)
names(bg.quad)[3] <- "monthYear"
bg.quad$dist <- 20
bg.quad$dist[bg.quad$site %in% "LUM3"] <- 10
bg.quad$ID <- createUniqueID(bg.quad, uniqueIDcols)



# site conditions (bay salinity, conductivity, temp)
# set bay conditions in distinct columns, rather than rows
siteCond <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_surfaceWater_151129.txt", skip = 2)
head(siteCond)
names(siteCond)[3] <- "monthYear"
siteCond$monthYear <- as.yearmon(siteCond$monthYear, "%b-%y")
siteCond$quadrat   <- "A"
siteCond$dist      <- 20
siteCond$dist[(siteCond$site %in% "LUM3")] <- 10
siteCond$ID        <- createUniqueID(siteCond, uniqueIDcols)

# split into bay, 5m and 20m
bay.cond    <- siteCond[siteCond$loc %in% c("BAY"), ]
five.cond   <- siteCond[siteCond$loc %in% c("5M"), ]
twenty.cond <- siteCond[siteCond$loc %in% c("20M"), ]

names(bay.cond)[4:6]    <- paste0(names(bay.cond)[4:6], ".bay")
names(five.cond)[4:8]   <- paste0(names(five.cond)[4:8], ".5m")
names(twenty.cond)[4:8] <- paste0(names(twenty.cond)[4:8], ".20m")




# soil properties
head(om2)
names(om2)[2] <- "monthYear"
om2$quadrat   <- "A"
om2$dist      <- 20
om2$dist[om2$site %in% "LUM3"] <- 10
om2$ID <- createUniqueID(om2, uniqueIDcols)


# benthic chl
chl2$quadrat     <- "A"
names(chl2)[1:2] <- c("site", "monthYear")
chl2$dist        <- 20
chl2$dist[chl2$site %in% "LUM3"] <- 10
chl2$ID          <- createUniqueID(chl2, uniqueIDcols)

# surface water nutrient data (porewater nutrients not available)
nuts.quad <- nuts[nuts$loc %in% c("BAY", "5M", "20M"), ]
names(nuts.quad)[8] <- "monthYear"
nuts.quad$quadrat   <- "A"
nuts.quad$dist      <- 20
nuts.quad$dist[(nuts.quad$site %in% "LUM3")] <- 10
nuts.quad$ID        <- createUniqueID(nuts.quad, uniqueIDcols)
nuts.quad[duplicated(nuts.quad$ID), ]

# split into bay, 5m and 20m
bay.nuts    <- nuts.quad[nuts.quad$loc %in% c("BAY"), ]
five.nuts   <- nuts.quad[nuts.quad$loc %in% c("5M"), ]
twenty.nuts <- nuts.quad[nuts.quad$loc %in% c("20M"), ]

names(bay.nuts)[2:5]    <- paste0(names(bay.nuts)[2:5], ".bay")
names(five.nuts)[2:5]   <- paste0(names(five.nuts)[2:5], ".5m")
names(twenty.nuts)[2:5] <- paste0(names(twenty.nuts)[2:5], ".20m")


gridC_plotData <- join_all(list(quadData, bg.quad[, c(5:7, 9)], 
                                om2[, c(3:13, 16)], chl2[, c(4:5, 18)], 
                                bay.cond[, c(4:6, 13)], bay.nuts[, c(1:5)], 
                                five.nuts[, 1:5], five.cond[, c(4:8, 13)],
                                twenty.nuts[, 1:5], twenty.cond[, c(4:8, 13)]),
                           by = "ID", type = "full")
gridC_plotData$site      <- sapply(strsplit(gridC_plotData$ID, split = "-"), "[", 1)
gridC_plotData$quadrat   <- sapply(strsplit(gridC_plotData$ID, split = "-"), "[", 2)
gridC_plotData$monthYear <- as.yearmon(sapply(strsplit(gridC_plotData$ID, split = "-"), "[", 3), "%b %Y")
gridC_plotData$month     <- substr(gridC_plotData$monthYear, 1, 3)
gridC_plotData$year      <- substr(gridC_plotData$monthYear, 5, 8)
gridC_plotData$dist      <- sapply(strsplit(gridC_plotData$ID, split = "-"), "[", 4)
for (i in 1:nrow(gridC_plotData)) {
  if (gridC_plotData[i, "site"] %in% paste0("LUM", 1:3)) {
    gridC_plotData$marsh[i] <- "LUMCON"    
  } else if (gridC_plotData[i, "site"] %in% c("TB1", "TB2")) {
    gridC_plotData$marsh[i] <- "Bay LaFleur"
  } else if (gridC_plotData[i, "site"] %in% c("TB3", "TB4")) {
    gridC_plotData$marsh[i] <- "Lake Barre"
  }
}
gridC_plotData <- gridC_plotData[order(gridC_plotData$marsh, gridC_plotData$site, gridC_plotData$monthYear, as.numeric(gridC_plotData$dist)), ]
gridC_plotData <- gridC_plotData[!is.na(gridC_plotData$marsh), ]
nrow(gridC_plotData)
length(unique(gridC_plotData$ID))

dataset12 <- gridC_plotData[, c("year", "month", "marsh", "latitude", "longitude", "site", "quadrat",
                                names(gridC_plotData)[6:12],
                                names(gridC_plotData)[17:19],
                                names(gridC_plotData)[31:32],
                                names(gridC_plotData)[c(20:22, 27:30, 36:39, 33:35)],
                                names(gridC_plotData)[c(40:57)]
                                )] # excludes porewater nutrients
# write.csv(dataset12, file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC1/plotData_160105.csv")

### aggregate data: annual data
abpp
# write.csv(abpp, file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC1/AnnualSiteData_151213.csv")

# or, very coarse: 
dd.napp # dd.napp.se
# write.csv(cbind(dd.napp, dd.napp.se[, -c(1:2)]), file = "C:/RDATA/SPAL_allometry/GRIIDC/CWC1/AnnualRegionData_151213.csv")


#####
##### 
#####
#####
