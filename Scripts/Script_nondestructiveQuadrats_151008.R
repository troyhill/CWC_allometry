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


### apply region-year-season-specific seasonal allometry (with year-specificity; i.e., "fall 2014 LUM" allometry would be used, rather than the average of all fall LUM1 parameters) to estimate stem masses
### region-season-specific parameters are in params[[1]]; plot-season-specific parameters are in plotParams[[1]]
hgts$mass <- as.numeric(NA)

for (i in 1:nrow(hgts)) {
  # identify site and season
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






#####
##### Or, apply plot-season-specific allometry
#####
# hgts$mass2 <- as.numeric(NA)
# 
# for (i in 1:nrow(hgts)) {
#   # identify plot and season
#   targetSite <- hgts$site[i]
#   targetSeas <- hgts$season[i]
#   
#   # find seasonal allometry parameters for the region
#   if (hgts$type[i] %in% "LIVE") {
#     coefficient <- plotParams[[1]]$coef.live[(plotParams[[1]]$site %in% targetSite) & (plotParams[[1]]$season %in% targetSeas)]
#     exponent    <- plotParams[[1]]$exp.live[(plotParams[[1]]$site %in% targetSite) & (plotParams[[1]]$season %in% targetSeas)]
#   } else if (hgts$type[i] %in% "DEAD") {
#     coefficient <- plotParams[[1]]$coef.dead[(plotParams[[1]]$site %in% targetSite) & (plotParams[[1]]$season %in% targetSeas)]
#     exponent    <- plotParams[[1]]$exp.dead[(plotParams[[1]]$site %in% targetSite) & (plotParams[[1]]$season %in% targetSeas)]
#   }
#   # apply allometry
#   hgts$mass2[i] <- coefficient * (hgts$hgt[i] ^ exponent) # bam!
# }
# 
# 
# summary(hgts$mass2)
# summary(hgts$mass2 - hgts$mass) # pretty minor difference, except 97 more NAs
# hist(hgts$mass2 - hgts$mass, 200, main = "", xlab = "Disparity between two mass estimates")
# 
# t.test(hgts$mass2, hgts$mass)
# 
# hgts[is.na(hgts$mass2), ]


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




