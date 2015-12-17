### this script estimates stem masses for the un-clipped quadrats using 
### the seasonal allometry equations from each region. Site-speciic allometry wasn't used because
### regional allometry is parameterized with larger sample sizes, and it may be the case that live/dead stems
### weren't found at a specific clip plot in a specific month, but were observed in the un-clipped quadrats
### order of business:
# 1: add season column to non-destructively measured plots
# 2: look up and apply allometry parameters

### This script is sourced in:
# Script_NAPPCalc_151006.R

### load litter data and get seasonal allometry parameters
source("C:/RDATA/SPAL_allometry/CWC_allometry/Scripts/Script_seasonalAllometry_151007.R")


### load stem heights
hgts <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_measuredStems_151020.txt")
head(hgts)

hgts <- marshName(hgts)
hgts       <- seasonLabel(hgts)
hgts$seas  <- substr(hgts$season, 1, 4)

# to calculate estimation errors
qc <- cwc

### Approach 1:
### apply region-year-season-specific allometry (pooling sites; i.e., "fall 2014 LUM" allometry would be used, rather than the average of all fall LUM1 parameters) to estimate stem masses
### region-season-specific parameters are in params[[1]]; plot-season-specific parameters are in plotParams[[1]]

for (i in 1:nrow(hgts)) {
  # identify region and season
  targetSite <- hgts$marsh[i]
  targetSeas <- hgts$season[i]
  
  # find seasonal allometry parameters for the region
    if (hgts$type[i] %in% "LIVE") {
      coefficient <- params[[1]]$coef.live[(params[[1]]$site %in% targetSite) & (params[[1]]$season %in% targetSeas)]
      exponent    <- params[[1]]$exp.live[(params[[1]]$site %in% targetSite) & (params[[1]]$season %in% targetSeas)]
    } else if (hgts$type[i] %in% "DEAD") {
      coefficient <- params[[1]]$coef.dead[(params[[1]]$site %in% targetSite) & (params[[1]]$season %in% targetSeas)]
      exponent    <- params[[1]]$exp.dead[(params[[1]]$site %in% targetSite) & (params[[1]]$season %in% targetSeas)]
    }
  # apply allometry
  hgts$mass[i] <- coefficient * (hgts$hgt[i] ^ exponent) # bam!
}

# for error calculation, model known masses

for (i in 1:nrow(qc)) {
  # identify region and season
  targetRegion <- qc$marsh[i]
  targetSeas <- qc$season[i]
  
  # find seasonal allometry parameters for the region
  if (qc$type[i] %in% "LIVE") {
    coefficient <- params[[1]]$coef.live[(params[[1]]$site %in% targetRegion) & (params[[1]]$season %in% targetSeas)]
    exponent    <- params[[1]]$exp.live[(params[[1]]$site %in% targetRegion) & (params[[1]]$season %in% targetSeas)]
  } else if (qc$type[i] %in% "DEAD") {
    coefficient <- params[[1]]$coef.dead[(params[[1]]$site %in% targetRegion) & (params[[1]]$season %in% targetSeas)]
    exponent    <- params[[1]]$exp.dead[(params[[1]]$site %in% targetRegion) & (params[[1]]$season %in% targetSeas)]
  }
  # apply allometry
  qc$mass_mod1[i] <- coefficient * (qc$hgt[i] ^ exponent) # bam!
}



##### Approach 2:
##### Apply region-season-specific allometry: seasonRegionParams
#####
seasonRegionParams$season <- substr(as.character(seasonRegionParams$time.period), nchar(as.character(seasonRegionParams$time.period)) - 3, nchar(as.character(seasonRegionParams$time.period)))
seasonRegionParams$region <- substr(as.character(seasonRegionParams$time.period), 1, nchar(as.character(seasonRegionParams$time.period)) - 5)
hgts$mass2 <- as.numeric(NA)

for (i in 1:nrow(hgts)) {
  # identify region and season
  targetRegion <- hgts$marsh[i]
  targetSeas <- substr(hgts$season[i], 1, 4)
  
  # find seasonal allometry parameters for the region
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- seasonRegionParams$coef.live[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
    exponent    <- seasonRegionParams$exp.live[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- seasonRegionParams$coef.dead[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
    exponent    <- seasonRegionParams$exp.dead[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
  }
  # apply allometry
  hgts$mass2[i] <- coefficient * (hgts$hgt[i] ^ exponent) # bam!
}

# calculate for known stem masses
for (i in 1:nrow(qc)) {
  # identify region and season
  targetRegion <- qc$marsh[i]
  targetSeas <- substr(qc$season[i], 1, 4)
  
  # find seasonal allometry parameters for the region
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- seasonRegionParams$coef.live[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
    exponent    <- seasonRegionParams$exp.live[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- seasonRegionParams$coef.dead[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
    exponent    <- seasonRegionParams$exp.dead[(seasonRegionParams$region %in% targetRegion) & (seasonRegionParams$season %in% targetSeas)]
  }
  # apply allometry
  qc$mass_mod2[i] <- coefficient * (qc$hgt[i] ^ exponent) # bam!
}





summary(hgts$mass2)
summary(hgts$mass2 - hgts$mass) # pretty minor difference, except 97 more NAs
hist(hgts$mass2 - hgts$mass, 200, main = "", xlab = "Disparity between two mass estimates")

t.test(hgts$mass2, hgts$mass)

hgts[is.na(hgts$mass2), ]



##### Approach 3:
##### Apply alternative region-season-specific allometry: 3 phenological seasons (season_altRegionParams)
#####
source("C:/RDATA/SPAL_allometry/Script_seasonality_explo_151205.R")
season_altRegionParams$season <- substr(as.character(season_altRegionParams$time.period), nchar(as.character(season_altRegionParams$time.period)) - 3, nchar(as.character(season_altRegionParams$time.period)))
season_altRegionParams$region <- substr(as.character(season_altRegionParams$time.period), 1, nchar(as.character(season_altRegionParams$time.period)) - 5)
hgts$mass3 <- as.numeric(NA)

for (i in 1:nrow(hgts)) {
  # identify region and season
  targetRegion <- hgts$marsh[i]
  if (substr(hgts$moYr[i], 1, 3) %in% seasons_alt[[1]]) {
    targetSeas <- names(seasons_alt)[1]
  } else if  (substr(hgts$moYr[i], 1, 3) %in% seasons_alt[[2]]) {
    targetSeas <- names(seasons_alt)[2]
  } else if  (substr(hgts$moYr[i], 1, 3) %in% seasons_alt[[3]]) {
    targetSeas <- names(seasons_alt)[3]
  }
  
  # find seasonal allometry parameters for the region
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- season_altRegionParams$coef.live[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
    exponent    <- season_altRegionParams$exp.live[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- season_altRegionParams$coef.dead[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
    exponent    <- season_altRegionParams$exp.dead[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
  }
  # apply allometry
  hgts$mass3[i] <- coefficient * (hgts$hgt[i] ^ exponent) # bam!
}


# calculate for known stem masses
for (i in 1:nrow(qc)) {
  # identify region and season
  targetRegion <- qc$marsh[i]
  if (substr(qc$time[i], 1, 3) %in% seasons_alt[[1]]) {
    targetSeas <- names(seasons_alt)[1]
  } else if  (substr(qc$time[i], 1, 3) %in% seasons_alt[[2]]) {
    targetSeas <- names(seasons_alt)[2]
  } else if  (substr(qc$time[i], 1, 3) %in% seasons_alt[[3]]) {
    targetSeas <- names(seasons_alt)[3]
  }
  
  # find seasonal allometry parameters for the region
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- season_altRegionParams$coef.live[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
    exponent    <- season_altRegionParams$exp.live[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- season_altRegionParams$coef.dead[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
    exponent    <- season_altRegionParams$exp.dead[(season_altRegionParams$region %in% targetRegion) & (season_altRegionParams$season %in% targetSeas)]
  }
  # apply allometry
  qc$mass_mod3[i] <- coefficient * (qc$hgt[i] ^ exponent) # bam!
}



summary(hgts$mass3)
summary(hgts$mass3 - hgts$mass)
hist(hgts$mass3 - hgts$mass, 200, main = "", xlab = "Disparity between two mass estimates")



##### Approach 4:
##### Apply single allometry equation (all data pooled)
#####
for (i in 1:nrow(hgts)) {
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- allParams$coef.live
    exponent    <- allParams$exp.live
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- allParams$coef.dead
    exponent    <- allParams$exp.dead
  }
  # apply allometry
  hgts$mass4[i] <- coefficient * (hgts$hgt[i] ^ exponent) # bam!
}

# calculate for known stem masses
for (i in 1:nrow(qc)) {
  # find seasonal allometry parameters for the region
  if (hgts$type[i] %in% "LIVE") {
    coefficient <- allParams$coef.live
    exponent    <- allParams$exp.live
  } else if (hgts$type[i] %in% "DEAD") {
    coefficient <- allParams$coef.dead
    exponent    <- allParams$exp.dead
  }
  # apply allometry
  qc$mass_mod4[i] <- coefficient * (qc$hgt[i] ^ exponent) # bam!
}



summary(qc$mass_mod1 - qc$mass)
summary(qc$mass_mod2 - qc$mass)
summary(qc$mass_mod3 - qc$mass) # not a dramatic effect, but of the 'seasonal' models this one performs better
summary(qc$mass_mod4 - qc$mass)

### calculate root mean squared error
rmse1 <- sqrt(sum((qc$mass_mod1 - qc$mass)^2, na.rm = T) / sum(!is.na(qc$mass)))
rmse2 <- sqrt(sum((qc$mass_mod2 - qc$mass)^2, na.rm = T) / sum(!is.na(qc$mass)))
rmse3 <- sqrt(sum((qc$mass_mod3 - qc$mass)^2, na.rm = T) / sum(!is.na(qc$mass)))
rmse4 <- sqrt(sum((qc$mass_mod4 - qc$mass)^2, na.rm = T) / sum(!is.na(qc$mass)))

#####
##### Figures
#####
# 
# ggplot(m.seaPlots[(m.seaPlots$variable %in% c("live.mean", "dead.mean")) & (m.seaPlots$marsh %in% ("LUM")), ], aes(x = seas, y = value)) + geom_point() + 
#   theme_bw() + theme(legend.title = element_blank()) + 
#   geom_errorbar(aes(ymin = value - se, ymax = value + se),
#                 width = 0) + 
# #   scale_colour_discrete(labels = c("LUMCON", "Bay La Fleur", "Lake Barre")) + 
#   labs(x = "Season", y = "Allometry exponent") +
#   facet_grid(. ~ variable, labeller = allomLabeller1)
# # ggsave("C:/RDATA/SPAL_allometry/LUM_allometryBySeason.png", width = 8, height= 6, units = "in", dpi = 300)




