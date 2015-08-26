library(plyr)
library(zoo)


proc_CWC_files <- function (dataset, 
                           ExcelParameterFile = "C:/RDATA/SPAL_allometry/data_LUM123/data_CWC_oldCoeffs_150821.txt") {
  # function returns parameters for exponential fits (and diagnostic plots, saved to working directory)
  # file structure changed in May 2015...
  if (names(dataset)[13] %in% "Total.Dried.Plant.Weight") {  
    names(dataset)[c(1:7, 13)] <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  } else names(dataset) <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  
  mo   <- as.character(grep(substr(dataset$time[1], 1, 3), month.abb))
  da   <- as.character(substr(dataset$time[1], 5, 6))
  moYr <- as.character(dataset$time[1])
  
  returnVals <- data.frame(month = as.character(dataset$time[1]), 
                           LUM1.R.coef = I(NA))
  
  
    
  for (i in 1:length(grep("LUM", levels(dataset$site)))) {
    x          <- dataset$hgt[(dataset$type %in% c("Live", "LIVE")) & (dataset$site %in% c(paste0("LUM", i), paste0("LUM ", i)))]
    y          <- dataset$mass[(dataset$type %in% c("Live", "LIVE")) & (dataset$site %in% c(paste0("LUM", i), paste0("LUM ", i)))]
    coefs      <- coef(model <- nls(y ~ I(a * exp(b * x)), start = list(a=0.1,b=0.1)))
    test       <- data.frame(xVals = 1:max(x))
    test$yVals <- coefs[1] * exp(coefs[2] * test$xVals)
    plotname   <- paste0("LUM", i, "-", moYr)

    png(filename = paste0(plotname, ".png"), width = 10, height = 8, units = "cm", res = 200)
      par(mar = c(4, 4, 0.3, 0.5), fig = c(0,1,0,1))
      plot(y ~ x, cex = 0.6, pch = 19, 
           xlab = "stem height (cm)",
           ylab = "stem mass (g)")
      text(x = as.numeric(quantile(na.omit(x), 0.4)), y = as.numeric(quantile(na.omit(y), 0.8)), labels = plotname)
      lines(y = test$yVals, x = test$xVals)
    dev.off()
    
    # calculate diagnostics
    predicted         <- coefs[1] * exp(coefs[2] * x)
    squared_error     <- (predicted - y)^2
    returnVals[paste0("LUM", i, ".R.coef")] <- coefs[1]
    returnVals[paste0("LUM", i, ".R.exp")]  <- coefs[2]
    returnVals[paste0("LUM", i, ".R.MSE")]  <- sum(squared_error, na.rm = T)
    returnVals[paste0("LUM", i, ".R.corr")] <- sqrt(1 - (deviance(model)  / sum((y[!is.na(y)] - mean(y, na.rm = T))^2)))
    
    
    oldParams <- read.delim(ExcelParameterFile, skip = 1)
    colNos <- list(c(2, 3), c(4, 5), c(6, 7))
    
    # find Excel coefficients by date matching
    returnVals[paste0("LUM", i, ".Excel.coef")] <- Excel.coef <- oldParams[oldParams$Month %in% moYr, colNos[[i]][1]]
    returnVals[paste0("LUM", i, ".Excel.exp")]  <- Excel.exp  <-oldParams[oldParams$Month %in% moYr, colNos[[i]][2]]
    
    # calculate Excel diagnostics
    E.predicted     <- Excel.coef * exp(Excel.exp * x)
    E.squared_error <- (E.predicted - y)^2
    returnVals[paste0("LUM", i, ".Excel.MSE")]  <- sum(E.squared_error,  na.rm = T)
    returnVals[paste0("LUM", i, ".Excel.corr")] <- cor(y[!is.na(y)], E.predicted[!is.na(y)])
    
    rm(coefs)
  }
  returnVals
}


getAllometryParams <- function (dataset, sitesIncluded = "all", 
                                returnData = "FALSE", combinePlots = "FALSE") {
  # function returns parameters for exponential fits (and diagnostic plots, saved to working directory)
  # This function is similar to proc_CWC_files, but differs in two important ways:
  # 1) it doesn't compare R and Excel, so only calculates R models. 
  # 2) This function also calculates allometric equations for dead material
  
  # dataset: a dataframe with plant masses, heights etc. 
  # sitesIncluded: which sites to report data for. "all" (default), "LUM", or "TB"
  # returnData: should raw data be returned? If "TRUE", 
  #  output is a list with two elements: allometry parameters and raw data.
  # combinePlots: if `TRUE`, data from all plots are pooled (if sitesIncluded = "all", all plots are pooled).
  
  # error checking
  countsAsTrue <- c("TRUE", "True", "T", "true")
  countsAsFalse <- c("FALSE", "False", "F", "false")
  TB.sites  <- paste0("TB", 1:4)
  LUM.sites <- paste0("LUM", 1:3)
  
  if (!sitesIncluded %in% c("all", "LUM", "TB")) {
    stop ("`sitesIncluded` argument is erroneous. Acceptable entries: `all`, `LUM`, `TB`")
  }
  if (length(sitesIncluded) != 1) {
    stop ("`sitesIncluded` argument should have one entry. Acceptable entries: `all`, `LUM`, `TB`")
  }
   if (!returnData %in% c(countsAsTrue, countsAsFalse)) {
    stop ("`returnData` argument is erroneous. Acceptable entries: `TRUE` or `FALSE`")
  }
  if (length(returnData) != 1) {
    stop ("`returnData` argument should have one entry. Acceptable entries: `TRUE` or `FALSE`")
  }
  if (!combinePlots %in% c(countsAsTrue, countsAsFalse)) {
    stop ("`combinePlots` argument is erroneous. Acceptable entries: `TRUE` or `FALSE`")
  }
  if (length(combinePlots) != 1) {
    stop ("`combinePlots` argument should have one entry. Acceptable entries: `TRUE` or `FALSE`")
  }
  
  # clean up the dataset
  if (names(dataset)[13] %in% "Total.Dried.Plant.Weight") {  
    names(dataset)[c(1:7, 13)] <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  } else names(dataset) <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
  
  # remove empty rows and columns
  dataset <- dataset[!is.na(dataset$ID), c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")]
  
  # sort dataframe to follow pattern LUM1-3, TB1-4
  dataset <- dataset[order(dataset$site), ]
  
  
  # homogenize inconsistent labeling
  dataset$site <- as.factor(gsub(pattern = " ", replacement = "", x = dataset$site))
  dataset$type <- as.character(dataset$type)
  
  # find month (robust to datasets spanning multiple months)
  mo <- paste(paste(as.character(grep(paste(c(unique(substr(dataset$time, 1, 3))), collapse = "|"), month.abb))), 
          collapse = " ")
  # year reports all years covered by a dataset; specific year-month combinations are reported by moYr
  da <- paste0(paste0("20", unique(as.character(substr(dataset$time, 5, 6)))), collapse = " ")
  moYr <- paste0(unique(as.character(dataset$time)), collapse = " ")
   

  # check site exclusion flag
  if ((sitesIncluded %in% "all") & (combinePlots %in% countsAsFalse)) {
    plotsInData <- levels(dataset$site)[levels(dataset$site) %in% c(TB.sites, LUM.sites)]
  } else if ((sitesIncluded %in% "all") & (combinePlots %in% countsAsTrue)) {
    # list all plots merged, if plots are being merged
    plotsInData <- paste0(levels(dataset$site)[levels(dataset$site) %in% c(TB.sites, LUM.sites)], collapse = " ")
  } else if ((sitesIncluded %in% "LUM") & (combinePlots %in% countsAsFalse)) {
    plotsInData <- levels(dataset$site)[levels(dataset$site) %in% LUM.sites]
    dataset <- dataset[dataset$site %in% LUM.sites, ]
  } else if ((sitesIncluded %in% "TB") & (combinePlots %in% countsAsFalse))  {
    plotsInData <- levels(dataset$site)[levels(dataset$site) %in% TB.sites]
    dataset <- dataset[dataset$site %in% TB.sites, ]
  } else if ((sitesIncluded %in% "LUM") & (combinePlots %in% countsAsTrue)) {
    plotsInData <- paste0(levels(dataset$site)[levels(dataset$site) %in% LUM.sites], collapse = " ")
    dataset <- dataset[dataset$site %in% LUM.sites, ]
  } else if ((sitesIncluded %in% "TB") & (combinePlots %in% countsAsTrue)) {
    plotsInData <- paste0(levels(dataset$site)[levels(dataset$site) %in% TB.sites], collapse = " ")
    dataset <- dataset[dataset$site %in% TB.sites, ]
  }
  
  if (length(plotsInData) == 0) {
    print(paste0("no plots from ", sitesIncluded, " were found on ", moYr, ". Try `all`."))
    returnVals <- data.frame(monthYear = moYr, 
                             month = mo, year = da,
                             site = NA, plot = NA,
                             coef.live = as.numeric(NA), exp.live = as.numeric(NA), 
                             MSE.live = as.numeric(NA), r.live = as.numeric(NA),    
                             coef.dead = as.numeric(NA), exp.dead = as.numeric(NA), 
                             MSE.dead = as.numeric(NA), r.dead = as.numeric(NA)
    )   
  } else {
    
  # determine how many plots are included
  siteOnly <- unlist(strsplit(as.character(dataset$site), split = c("[1-9]")))
  sitesInData <- unique(siteOnly)

  # set number of rows in allometry parameter file
  # if all sites are selected, this step begins to combine  their data
  if ((combinePlots %in% countsAsTrue) & (sitesIncluded %in% "all")) {
    n <- 1
  } else if ((combinePlots %in% countsAsTrue) & (!sitesIncluded %in% "all")) {
    n <- 1
  } else {
    n <- length(plotsInData) 
  }
  
  returnVals <- data.frame(monthYear = rep(moYr, times = n), 
                           month = rep(mo, times = n), year = rep(da, times = n),
                           site = rep(NA, times = n), plot = rep(NA, times = n),
                           coef.live = rep(as.numeric(NA), times = n), exp.live = rep(as.numeric(NA), times = n), # y = coef * e^(exp * b)
                           MSE.live = rep(as.numeric(NA), times = n), r.live = rep(as.numeric(NA), times = n),    # diagnostics
                           coef.dead = rep(as.numeric(NA), times = n), exp.dead = rep(as.numeric(NA), times = n), 
                           MSE.dead = rep(as.numeric(NA), times = n), r.dead = rep(as.numeric(NA), times = n)
                           )    
  
  for (i in 1:n) {
    
    if (combinePlots %in% countsAsTrue) {
      siteName <- strsplit(plotsInData, " ")[[1]]
    }  else {
      siteName <- plotsInData[i]
    }
    
    if ((length(sitesInData) > 1) & (combinePlots %in% countsAsTrue)) {
      returnVals$site[i] <- paste0(sitesInData, collapse = " ")
    } else {
      returnVals$site[i] <- strsplit(as.character(unique(dataset$site)), split = c("[1-9]| "))[[i]][1] # "LUM" or "TB"
    }
    
    # `plot` reports `sitesIncluded` if data are being combined across plots
    if (combinePlots %in% countsAsFalse) {
      returnVals$plot[i]  <- plotsInData[i]
    } else if ((combinePlots %in% countsAsTrue) & (sitesIncluded %in% "all")) {
      returnVals$plot[i]  <- plotsInData[i]
    } else if ((combinePlots %in% countsAsTrue) & (!sitesIncluded %in% "all")) {
      returnVals$plot[i] <-  returnVals$site[i]
    }
    
    ### data subset by type
    x          <- dataset$hgt[(dataset$type %in% c("Live", "LIVE", "live")) & (dataset$site %in% siteName)]
    y          <- dataset$mass[(dataset$type %in% c("Live", "LIVE", "live")) & (dataset$site %in% siteName)]
    x.dead     <- dataset$hgt[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]
    y.dead     <- dataset$mass[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]
    
### This didn't work, not sure why
#     if one x-vector has length zero but the other doesn't, it means the plot was sampled and 
#     both types should be in the data (if data is returned by function)
#     if (returnData %in% countsAsTrue) {
#       # if one type has data but the other doesn't insert 0s in dataset (because plot was sampled)
#       if ((length(x) == 0) & (length(x.dead) != 0)) { 
#         fillData <- dataset[1, ]
#         fillData$type <- "LIVE"
#         fillData$ID <- fillData$tin <- fillData$tin_plant <- NA
#         fillData$hgt <- fillData$mass <- as.numeric(0)
#         dataset <- rbind(dataset, fillData)
#       } else if ((length(x.dead) == 0) & (length(x) != 0)) {
#         fillData <- dataset[1, ]
#         fillData$type <- "DEA"
#         fillData$ID <- fillData$tin <- fillData$tin_plant <- NA
#         fillData$hgt <- fillData$mass <- as.numeric(0)
#         dataset <- rbind(dataset, fillData)
#         }
#     }
###
    
    
    if (length(x) < 2) {
      returnVals$coef.live[i] <- as.numeric(NA)
      returnVals$exp.live[i]  <- as.numeric(NA)
      returnVals$MSE.live[i]  <- as.numeric(NA)
      returnVals$r.live[i]    <- as.numeric(NA)
    } else {
      coefs      <- coef(model <- nls(y ~ I(a * exp(b * x)), start = list(a = 0.1, b = 0.1)))
      predicted           <- coefs[1] * exp(coefs[2] * x)
      squared_error       <- (predicted - y)^2
      returnVals$coef.live[i] <- coefs[1]
      returnVals$exp.live[i]  <- coefs[2]
      returnVals$MSE.live[i]  <- sum(squared_error, na.rm = T)
      returnVals$r.live[i]    <- sqrt(1 - (deviance(model)  / sum((y[!is.na(y)] - mean(y, na.rm = T))^2)))
    }
    
    # do the same for dead stems, if there's more than two stems
    if (length(dataset$hgt[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]) > 2) {
        if (length(x.dead) < 2) {
          returnVals$coef.dead[i] <- as.numeric(NA)
          returnVals$exp.dead[i]  <- as.numeric(NA)
          returnVals$MSE.dead[i]  <- as.numeric(NA)
          returnVals$r.dead[i]    <- as.numeric(NA)
        } else {
        coefs.dead <- coef(model.dead <- nls(y.dead ~ I(a * exp(b * x.dead)), start = list(a = 0.1, b = 0.1)))
        
        predicted.dead      <- coefs.dead[1] * exp(coefs.dead[2] * x.dead)
        squared_error.dead  <- (predicted.dead - y.dead)^2
        
        
        returnVals$coef.dead[i] <- coefs.dead[1]
        returnVals$exp.dead[i]  <- coefs.dead[2]
        returnVals$MSE.dead[i]  <- sum(squared_error.dead, na.rm = T)
        returnVals$r.dead[i]    <- sqrt(1 - (deviance(model.dead)  / sum((y.dead[!is.na(y.dead)] - mean(y.dead, na.rm = T))^2)))
        } 
    }
  }
  }
  if (returnData %in% countsAsFalse) {
    returnVals
  } else if (returnData %in% countsAsTrue) {
    returnVals <- list(parameters = returnVals, data = dataset)
  } 
}


batch <- function (inputList, fun, ...) {
  # function allows batch analysis of a list of dataframes. Output is a dataframe of results
  fun <- match.fun(fun)
  for(j in 1:length(inputList)) {
    a <- fun(inputList[[j]], ...)
    if (j != 1) {
      if (is.list(a)) {
        outputObj <- mapply(rbind, outputObj, a, SIMPLIFY = FALSE)
      } else {
        outputObj <- rbind(outputObj, a)
      }
    } else {
      outputObj <- a
    }
  }
  if (is.list(a)) {
    rownames(outputObj[[2]]) <- 1:nrow(outputObj[[2]])
  }
  outputObj
}




plotAllom <- function(monthlyData, site = "LUM", type = "both", save = "TRUE") {
  procData <- getAllometryParams(dataset = monthlyData, returnData = "TRUE", sitesIncluded = site)
  
  living <- c("Live", "coef.live", "exp.live", "live")
  dying <- c("Dead", "coef.dead", "exp.dead", "dead")

  if (save %in% c("TRUE", "True", "true", "T")) {
    png(filename = paste0("Allom-", procData[[2]]$time, ".png"), width = 15, height = 8, units = "cm", res = 300)
  }
  if (type %in% c("live", "Live", "LIVE", "both", "Both", "BOTH")) {
     if (type %in% c("both", "Both", "BOTH"))  {
      par(mar = c(4, 4, 0.3, 0.5), fig = c(0, 0.48, 0, 1))
    } else if (!type %in% c("both", "Both", "BOTH")) {
      par(mar = c(4, 4, 0.3, 0.5), fig = c(0, 1, 0, 1))
    }
    plot(procData[[2]]$mass ~ procData[[2]]$hgt,
         ylab = "mass (g)", xlab = "height (cm)",
         type = "n", las = 1)
    
    for (i in 1:nrow(procData[[1]])) {
      points(x = procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% living[1])], 
             y = procData[[2]]$mass[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% living[1])], 
             pch = i, cex = 0.8, las = 1)
      # predicted values
      min.x <- min(procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% living[1])], na.rm = T)
      max.x <- max(procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% living[1])], na.rm = T)
      xVals <- c((min.x * 100):(max.x * 100)) / 100
      modeled <- procData[[1]][paste0("coef.", living[4])][[1]][i] * exp(procData[[1]][paste0("exp.", living[4])][[1]][i] *
                                                    xVals)
      lines(x = xVals, y = modeled, lty = i)
    }
    legend(x = 1.1 * min(procData[[2]]$hgt, na.rm = T), y = 0.7 * max(procData[[2]]$mass, na.rm = T), 
           legend = c(unique(procData[[1]]$plot)), pch = c(1:nrow(procData[[1]])), 
           lty = c(1:nrow(procData[[1]])),
           cex = 0.7,  merge = TRUE, bty = "n", 
           title = paste0(procData[[1]]$monthYear[1], " ", living[4], " stems")
    )
  }
  if (type %in% c("dead", "Dead", "DEAD", "both", "Both", "BOTH")) {
    if (type %in% c("both", "Both", "BOTH"))  {
      par(new = TRUE)
      par(mar = c(4, 4, 0.3, 0.5), fig = c(0.5, 1, 0, 1))
    } else if (!type %in% c("both", "Both", "BOTH")) {
      par(mar = c(4, 4, 0.3, 0.5), fig = c(0, 1, 0, 1))
    }
    plot(procData[[2]]$mass ~ procData[[2]]$hgt,
         ylab = "mass (g)", xlab = "height (cm)",
         type = "n", las = 1)
    
    for (i in 1:nrow(procData[[1]])) {
      points(x = procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% dying[1])], 
             y = procData[[2]]$mass[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% dying[1])], 
             pch = i, cex = 0.8, las = 1)
      # predicted values
      min.x <- min(procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% dying[1])], na.rm = T)
      max.x <- max(procData[[2]]$hgt[(procData[[2]]$site %in% procData[[1]]$plot[i]) & (procData[[2]]$type %in% dying[1])], na.rm = T)
      xVals <- c((min.x * 100):(max.x * 100)) / 100
      modeled <- procData[[1]][paste0("coef.", dying[4])][[1]][i] * 
                  exp(procData[[1]][paste0("exp.", dying[4])][[1]][i] * xVals)
      lines(x = xVals, y = modeled, lty = i)
    }
    legend(x = 1.1 * min(procData[[2]]$hgt, na.rm = T), y = 0.7 * max(procData[[2]]$mass, na.rm = T), 
           legend = c(unique(procData[[1]]$plot)), pch = c(1:nrow(procData[[1]])), 
           lty = c(1:nrow(procData[[1]])),
           cex = 0.7,  merge = TRUE, bty = "n", 
           title = paste0(procData[[1]]$monthYear[1], " ", dying[4], " stems")
    )
  }
  if (save %in% c("TRUE", "True", "true", "T")) {
    dev.off()
  }
}




PSC <- function(dataset, liveCol = "live", deadCol = "dead", yearCol = "year", siteCol = "site", type = "both") {
  # peak standing crop estimates of NAPP
  # runs for entire dataset, reports results by year for each site (summary statistics)
  #   dataset = dataframe with your data
  #   liveCol = name of the column with live biomass data
  #   deadCol = name of the column with dead biomass data
  #   yearCol = name of the column with year (4 digits)
  #   siteCol = name of the column with site name
  #   type = "PSC-A" is the traditional approach, using just live material; "PSC-B" uses live + dead material in estimates of peak standing crop
  
  # Usage examples: 
  # A single site, single year
  #   test <- smalley.prep[(smalley.prep$site %in% "LUM1") & (smalley.prep$year %in% "2014"), 1:6]
  #   PSC(test)
  # Single site, multiple years
  #   test2 <- smalley.prep[(smalley.prep$site %in% "LUM1"), 1:6]
  #   PSC(test2)
  # Multiple sites, multiple years
  #   test3 <- smalley.prep[, 1:6]
  #   PSC(test3)
  
  ### error checking
  if (sum(c(liveCol, deadCol, yearCol, siteCol) %in% names(dataset)) < 4) {
    stop ("Check column names. One or more column names were not found in the dataset.")
  }
  if (!type %in% c("both", "PSC-A", "PSC-B")) {
    stop ("`type` argument must be either `PSC-A`, `PSC-B`, or `both`")
  }
  if (length(type) > 1) {
    stop ("`type` argument is too long. It must be of length 1.")
  }
  ###
  
  for (h in 1:length(unique(dataset[, siteCol]))) {
    targetSite <- unique(dataset[, siteCol])[h]
    subData1 <- dataset[dataset[, siteCol] %in% targetSite, ]
    for (i in 1:length(unique(subData1[, yearCol]))) {
      targetYear <- unique(subData1[, yearCol])[i]
      subData <- subData1[subData1[, yearCol] %in% targetYear, ]
      if (type %in% "both") {
        PSC_A  <- max(subData[, liveCol], na.rm = T)
        PSC_B  <- max((subData[, liveCol] + subData[, deadCol]), na.rm = T)
        output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A, psc_b = PSC_B)
      } else if (type %in% "PSC-A") {
        PSC_A <- max(subData[, liveCol], na.rm = T)
        output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A)
      } else if (type %in% "PSC-B") {
        PSC_B  <- max((subData[, liveCol] + subData[, deadCol]), na.rm = T)
        output <- data.frame(site = targetSite, year = targetYear, psc_b = PSC_B)
      }
      
      if (i == 1) {
        intOutput <- output
      } else {
        intOutput <- rbind(intOutput, output)
      }
    }
    if (h == 1) {
      finalOutput <- intOutput
    } else {
      finalOutput <- rbind(finalOutput, intOutput)
    }
  }
  finalOutput
}



nappCalc <- function(dataset, liveCol = "live", deadCol = "dead", yearCol = "year", siteCol = "site", 
                          MH = "TRUE", summarize = "FALSE", timeCol = "time") {
  # implements Smalley (1959) and Millner and Hughes (1968)
  # runs for entire dataset, reports results by year for each site
  #   dataset = dataframe with your data
  #   liveCol = name of the column with live biomass data
  #   deadCol = name of the column with dead biomass data
  #   yearCol = name of the column with year (4 digits)
  #   siteCol = name of the column with site name
  #   timeCol = column with time data (in form "%B %Y"). This only matters for the summary statistics, which report peak timing
  #   MH      = if "TRUE", also implements Millner & Hughes 1968 (sum of positive changes in standing live biomass)
  #   summarize = if "TRUE", summary statistics (max NAPP estimates) are reported. TODO: report peak timing
  #
  # Usage examples: 
  # # Single site, single year
  #   test <- smalley.prep[(smalley.prep$site %in% "LUM1") & (smalley.prep$year %in% "2014"), 1:6]
  #   nappCalc(test)
  # # Single site, multiple years
  #   test2 <- smalley.prep[(smalley.prep$site %in% "LUM1"), 1:6]
  #   nappCalc(test2)
  # # Multiple sites, multiple years
  #   test3 <- smalley.prep[, 1:6]
  #   nappCalc(test3)
  #   nappCalc(test3, summarize = "TRUE")[[2]]
  # # Last two columns are identical to 
  #   PSC(smalley.prep[, 1:6])
  
  ### error checking
  countsAsTrue  <- c("T", "TRUE", "true", "True")
  countsAsFalse <- c("F", "FALSE", "false", "False")
  
  if (!MH %in% c(countsAsTrue, countsAsFalse)) {
    stop ("`MH` argument isn't recognized. Input can be either `TRUE` or `FALSE`.")
  }
  if (sum(c(liveCol, deadCol, yearCol, siteCol) %in% names(dataset)) < 4) {
    stop ("Check column names. One or more column names were not found in the dataset.")
  }
  ###
  
  tempData <- dataset
  # column names as variables:
  smalley.inc <- "smalley.inc"
  smalley     <- "smalley"
  live.inc    <- "live.inc" # VTS1975's delL
  dead.inc    <- "dead.inc" # VTS1975's delD
  MH          <- "MillnerHughes"
  # Valiela, Teal, Sass 1975
  eV          <- "VTS1975.inc"
  VTS         <- "VTS1975"
  PSC_A       <- "psc.live"
  PSC_B       <- "psc.tot"
  
  
  # more variables than necessary are appended to dataset
  tempData[, PSC_A] <- tempData[, PSC_B] <- tempData[, VTS] <- tempData[, MH] <- tempData[, smalley] <- tempData[, smalley.inc] <- 
    tempData[, eV] <- tempData[, dead.inc] <- tempData[, live.inc] <- as.numeric(NA) 
  
  for (h in 1:length(unique(tempData[, siteCol]))) {
    targetSite <- unique(tempData[, siteCol])[h]
    subData1 <- tempData[tempData[, siteCol] %in% targetSite, ]
    # calculates biomass increments (Bn - B(n-1)) over all available years
    for (j in 2:(nrow(subData1))) {
      # live biomass increments
      subData1[, live.inc][j] <- subData1[, liveCol][j] - subData1[, liveCol][j - 1]
      # dead biomass increments
      subData1[, dead.inc][j] <- subData1[, deadCol][j] - subData1[, deadCol][j - 1]
      
      # calculate e from Valiela, Teal, Sass 1975
      if ((subData1[, live.inc][j] >= 0) & (subData1[, dead.inc][j] < 0)) {
        subData1[, eV][j] <- -subData1[, dead.inc][j]
      } else if (subData1[, live.inc][j] < 0) {
        subData1[, eV][j] <- -(subData1[, dead.inc][j] + subData1[, live.inc][j])
      }
      # change e to zero, if negative
      if (!is.na(subData1[, eV][j])) {
        if (subData1[, eV][j] < 0 ) {
          subData1[, eV][j] <- 0
        }
      }
      
      # apply decision rules
        # both increments positive: sum of both
        # both negative: zero
        # live + and dead -: use live
        # live - and dead +: difference (if positive, otherwise use zero)
      if ((subData1[, live.inc][j] > 0) & (subData1[, dead.inc][j] > 0)) {
        newVal <- subData1[, live.inc][j] + subData1[, dead.inc][j]
      } else if ((subData1[, live.inc][j] <= 0) & (subData1[, dead.inc][j] <= 0)) {
        newVal <- 0
      } else if ((subData1[, live.inc][j] > 0) & (subData1[, dead.inc][j] <= 0)) {
        newVal <- subData1[, live.inc][j] 
      } else if ((subData1[, live.inc][j] <= 0) & (subData1[, dead.inc][j] > 0)) {
        calc <- subData1[, live.inc][j] + subData1[, dead.inc][j]
        if (calc < 0) {
          newVal <- 0
        } else {
          newVal <- calc
        }
      }
      # write values to output dataframe
      tempData[tempData[, siteCol] %in% targetSite, live.inc][j]    <- subData1[, live.inc][j]
      tempData[tempData[, siteCol] %in% targetSite, dead.inc][j]    <- subData1[, dead.inc][j]
      tempData[tempData[, siteCol] %in% targetSite, eV][j]          <- subData1[, eV][j]
      tempData[tempData[, siteCol] %in% targetSite, smalley.inc][j] <- newVal
    }
      
  # now, calculate NAPP cumulatively for each year
  # treat NAs as zeroes in NAPP calculation
  for (k in 1:length(unique(subData1[, yearCol]))) {
    targetYear <- unique(subData1[, yearCol])[k]
    # data for a single site, single year
    subData2   <- tempData[(tempData[, yearCol] %in% targetYear) & (tempData[, siteCol] %in% targetSite), ]
    
    # sum smalley increments
    subData2[, smalley][!is.na(subData2[, smalley.inc])] <- cumsum(subData2[, smalley.inc][!is.na(subData2[, smalley.inc])])
    
    # sum live biomass increments for Millner & Hughes 1968
    temp <- subData2[, live.inc][subData2[, live.inc] > 0][!is.na(subData2[, live.inc][subData2[, live.inc] > 0])]
    # this line could be problematic if two live biomass increments are identical in a single year
    subData2[which(subData2[, live.inc] %in% temp), MH] <- cumsum(temp[!is.na(temp)])
    
    # sum e from Valiela, Teal, Sass 1975
    subData2[, VTS][!is.na(subData2[, eV])] <- cumsum(subData2[, eV][!is.na(subData2[, eV])])
      
    # sum peak standing crop increments
    subData2[, PSC_A][!is.na(subData2[, liveCol])]                       <- cummax(subData2[, liveCol][!is.na(subData2[, liveCol])])
    subData2[, PSC_B][!is.na(subData2[, liveCol] + subData2[, deadCol])] <- cummax(c(subData2[, liveCol] + subData2[, deadCol])[!is.na(c(subData2[, liveCol] + subData2[, deadCol]))])
    
    
    # add to output dataframe
    tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), smalley]  <- subData2[, smalley]
    tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), MH]       <- subData2[, MH]
    tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), VTS]      <- subData2[, VTS]
    tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_A]    <- subData2[, PSC_A]
    tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_B]    <- subData2[, PSC_B]
    }
  }
  
  if (summarize %in% countsAsFalse) {
    output <- tempData
  } else if (summarize %in% countsAsTrue) {
    # this is a weak point; I'm not positive this will work if column names differ
    summaryStats <- ddply(tempData, .(eval(parse(text = siteCol)), eval(parse(text = yearCol))), summarise,
                          napp.smalley   = max(eval(parse(text = smalley)), na.rm = T),
                          napp.MH        = max(eval(parse(text = MH)), na.rm = T),
                          napp.VTS       = max(eval(parse(text = VTS)), na.rm = T),
                          napp.psc.a     = max(eval(parse(text = PSC_A)), na.rm = T),
                          napp.psc.b     = max(eval(parse(text = PSC_B)), na.rm = T)
                          
                          # Couldn't get the peak timing for smalley to work...
                          #  ddply(tempData, .(eval(parse(text = siteCol)), eval(parse(text = yearCol))), summarise,
                          #         napp.smalley   = max(eval(parse(text = smalley)), na.rm = T),
                          #         napp.MH        = max(eval(parse(text = MH)), na.rm = T),
                          #         timing.MH      = na.omit(eval(parse(text = timeCol))[eval(parse(text = MH))      == max(eval(parse(text = MH)), na.rm = T)])[1],
                          #         timing.smalley = na.omit(eval(parse(text = timeCol))[eval(parse(text = smalley)) == napp.smalley])[1]
                          #         )
                          # na.omit(tempData[tempData[, MH] == max(tempData[, MH], na.rm = T), timeCol])[1]
                          # na.omit(tempData[tempData[, smalley] == max(tempData[, smalley], na.rm = T), timeCol])[1] # works
    )
    names(summaryStats)[1:2] <- c( "site", "year")
    
#     # get timing manually
#     # includes multiple times per year
#     timing.MH      <- data.frame(tempData[, timeCol][which(tempData[, MH] %in% summaryStats$napp.MH)])
#     timing.MH$year <- substr(timing.MH[, 1], 5, 9)
#     summaryStats$timing.MH <- summaryStats$timing.MH <- NA
#     for(i in 1:length(unique(timing.MH$year))) {
#       # use first peak
#       peakTime   <- as.character(timing.MH[timing.MH$year %in% unique(timing.MH$year)[i], 1][1]) # could switch to as.yearmon
#       insertYear <- timing.MH[timing.MH$year %in% unique(timing.MH$year)[i], "year"][1]
#       summaryStats$timing.MH[summaryStats$year %in% insertYear] <- peakTime
#     }
#     
#     timing.smalley      <- data.frame(tempData[, timeCol][which(tempData[, smalley] %in% summaryStats$napp.smalley)])
#     timing.smalley$year <- substr(timing.smalley[, 1], 5, 9)
#     for(i in 1:length(unique(timing.smalley$year))) {
#       peakTime   <- as.character(timing.smalley[timing.smalley$year %in% unique(timing.smalley$year)[i], 1][1]) # could switch to as.yearmon
#       insertYear <- timing.smalley[timing.smalley$year %in% unique(timing.smalley$year)[i], "year"][1]
#       summaryStats$timing.smalley[summaryStats$year %in% insertYear] <- peakTime
#     }
    output <- list(intervalData = tempData, summary = summaryStats)
  }
  output
}





