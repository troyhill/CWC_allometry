CWC_aggregate <- function(dataset) {
  # function aggregates CWC data to estimate per-plot biomass, stem density, and average longest-culm length
  TB.sites  <- paste0("TB", 1:4)
  LUM.sites <- paste0("LUM", 1:3)
  plotSize <- 0.25^2 # square meters
  
  if (names(dataset)[13] %in% "Total.Dried.Plant.Weight") {  
    names(dataset)[c(1:7, 13)] <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  } else names(dataset) <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  
  # remove empty rows and columns
  dataset <- dataset[!is.na(dataset$ID), c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")]
  
  # sort dataframe to follow pattern LUM1-3, TB1-4
  dataset <- dataset[order(dataset$site), ]
  
  
  # homogenize inconsistent labeling
  dataset$site <- as.factor(gsub(pattern = " ", replacement = "", x = dataset$site))
  
  # find month (robust to datasets spanning multiple months)
  mo <- paste(paste(as.character(grep(paste(c(unique(substr(dataset$time, 1, 3))), collapse = "|"), month.abb))), 
              collapse = " ")
  # year reports all years covered by a dataset; specific year-month combinations are reported by moYr
  da <- paste0(paste0("20", unique(as.character(substr(dataset$time, 5, 6)))), collapse = " ")
  moYr <- paste0(unique(as.character(dataset$time)), collapse = " ")
  
  # determine how many plots are included
  siteOnly <- unlist(strsplit(as.character(dataset$site), split = c("[1-9]")))
  sitesInData <- levels(as.factor(siteOnly)) 
  
  # set number of rows in allometry parameter file
  # if all sites are selected, this step begins to combine  their data
  plotsInData <- levels(dataset$site)[levels(dataset$site) %in% c(TB.sites, LUM.sites)]
  n <- length(plotsInData)
  
  returnVals <- data.frame(year   = rep(da, times = n), 
                           month   = rep(mo, times = n), 
                           monthYear = rep(moYr, times = n), 
                           plot = rep(as.numeric(NA), times = n),        biomass.g.m2 = rep(as.numeric(NA), times = n),
                           dead.biomass.g.m2 = rep(as.numeric(NA), times = n), 
                           stemDensity = rep(as.numeric(NA), times = n), dead.stemDensity = rep(as.numeric(NA), times = n), 
                           longestThree = rep(as.numeric(NA), times = n))
  
  for (i in 1:n) {
    returnVals$plot[i]              <- siteName <- plotsInData[i]
    newData                         <- dataset[(dataset$site %in% siteName), c("type", "hgt", "mass")]
#     newData                         <- newData[rowSums(is.na(newData)) != ncol(newData), ] # remove rows with all NAs
    returnVals$biomass.g.m2[i]      <- sum(newData$mass[newData$type %in% c("Live", "LIVE", "live")], na.rm = T)  / plotSize
    returnVals$dead.biomass.g.m2[i] <- sum(newData$mass[newData$type %in% c("Dead", "DEAD", "dead")], na.rm = T)  / plotSize
    returnVals$stemDensity[i]       <- nrow(newData[as.character(newData$type) %in% c("Live", "LIVE", "live"), ]) / plotSize # live stems
    returnVals$dead.stemDensity[i]  <- nrow(newData[as.character(newData$type) %in% c("Dead", "DEAD", "dead"), ]) / plotSize # dead stems
    returnVals$longestThree[i]      <- mean(sort(newData$hgt[newData$type %in% c("Live", "LIVE", "live")], decreasing = TRUE)[1:3], na.rm = T)  
  }
  returnVals[!is.na(returnVals$plot), ]
}

