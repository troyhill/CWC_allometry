library(plyr)
library(zoo)

marshName <- function(data, siteCol = "site") {
  # function takes a dataset and the name of the column with site names (e.g., "LUM1"), and 
  # adds a column with marsh names ("LUM", "TB-A", "TB-B")
  tempData <- data
  tempData$marsh <- as.character(NA)
  for (i in 1:nrow(tempData)) {
    if (tempData[i, siteCol] %in% paste0("LUM", 1:3)) {
      tempData$marsh[i] <- "LUM"    
    } else if (tempData[i, siteCol] %in% c("TB1", "TB2")) {
      tempData$marsh[i] <- "TB-A"
    } else if (tempData[i, siteCol] %in% c("TB3", "TB4")) {
      tempData$marsh[i] <- "TB-B"
    }
  }
  tempData
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
  } else {
    rownames(outputObj)      <- 1:nrow(outputObj)
  }
  outputObj
}

mergeMonths <- function (dataset_list) {
  # function formats and merges monthly tab-delimited clip-plot data files from 2014 and 2015
  # dataset_list: a dataframe with plant masses, heights etc. 
  
  for(j in 1:length(dataset_list)) {
    dataset <- dataset_list[[j]]
    # clean up the dataset
    if (names(dataset)[13] %in% "Total.Dried.Plant.Weight") {  
      names(dataset) <- c("site", "time", "type", "ID", "hgt", "tin", "stm_ind", "stemMass", 
                          "lf_ind", "leafMass", "thirdLeaf_ind", "thirdLeafMass",
                          "mass")
      # re-calculate leaf/stem masses, subtracting tin mass
      dataset$leafMass      <- dataset$leafMass - (dataset$lf_ind * dataset$tin)
      dataset$thirdLeafMass <- dataset$thirdLeafMass - (dataset$thirdLeaf_ind * dataset$tin)
      dataset$stemMass      <- dataset$stemMass - (dataset$stm_ind * dataset$tin)
    } else if (!names(dataset)[8] %in% "mass") {
      names(dataset) <- c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")
      dataset[, c("stm_ind", "stemMass", "lf_ind", "leafMass", "thirdLeaf_ind", "thirdLeafMass")] <- as.numeric(NA)
    }
    
    # "leafMass" column has combined leaf masses, where available.
    dataset$leafMass <- ifelse(!is.na(dataset$thirdLeafMass), dataset$leafMass + dataset$thirdLeafMass, dataset$leafMass)
    
    # remove empty rows and columns (determined by empty tin ID)
    dataset          <- dataset[!is.na(dataset$ID), c("site", "time", "type", "ID", "hgt", 
                                                      "stemMass", "leafMass", "mass")]
    
    # sort dataframe to follow pattern LUM1-3, TB1-4
    dataset          <- dataset[order(dataset$site), ]
    
    # homogenize inconsistent labeling
    dataset$site <- as.factor(gsub(pattern = " ", replacement = "", x = as.character(dataset$site)))
    dataset$type <- as.character(dataset$type)
    dataset$ID   <- as.character(dataset$ID)
    dataset$time <- as.character(dataset$time)
    
    if (j != 1) {
      outputDF <- rbind(outputDF, dataset)
    } else {
      outputDF <- dataset
    }
  }
  outputDF
}

proc_CWC_files <- function (dataset, startVals = 0.2,
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
    coefs      <- coef(model <- nls(y ~ I(a * x^b), start = list(a = startVals, b = startVals)))
    test       <- data.frame(xVals = 1:max(x))
    test$yVals <- coefs[1] * test$xVals^coefs[2]
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
    returnVals[paste0("LUM", i, ".R.corr")] <- 1 - (deviance(model)  / sum((y[!is.na(y)] - mean(y, na.rm = T))^2))
    
    
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


getParams <- function (dataset, massCol = "mass", heightCol = "hgt", typeCol = "type", 
                       cutOff = 5, startVals = 0.2,
                       timeName = NA, dataName = NA,
                       returnPlot = FALSE, savePlot = FALSE, plotTitle = NA) {
  # returnPlot: live stems are on the left, dead are on right
  # typeCol: column indicating live/dead status
  # a stream-lined workhorse of getAllometryParams. Generates live/dead allometry parameters for whatever dataset it's fed. 
  # really, this is just a wrapper for nls().
  # getParams(cwc[(cwc$season %in% "sprg 13") & (cwc$marsh %in% "LUM"), ])
  # getParams(cwc[(cwc$season %in% "sprg 13") & (cwc$marsh %in% "LUM"), ], returnPlot = TRUE)
  # getParams(cwc[(cwc$season %in% "sprg 13") & (cwc$marsh %in% "LUM"), ], returnPlot = TRUE, plotTitle = "sprg 13")
  
  returnVals <- data.frame(time.period = timeName, 
                           data.name = dataName,
                           coef.live = as.numeric(NA), exp.live = as.numeric(NA),
                           MSE.live = as.numeric(NA), r.live = as.numeric(NA),
                           coef.dead = as.numeric(NA), exp.dead = as.numeric(NA),
                           MSE.dead = as.numeric(NA), r.dead = as.numeric(NA)
  )
  
  x          <- dataset[, heightCol][tolower(dataset[, typeCol]) %in% "live"]
  y          <- dataset[, massCol][tolower(dataset[, typeCol]) %in% "live"]
  x.dead     <- dataset[, heightCol][tolower(dataset[, typeCol]) %in% "dead"]
  y.dead     <- dataset[, massCol][tolower(dataset[, typeCol]) %in% "dead"]
  
  # get live biomass parameters
  if (length(x) < cutOff) {
    returnVals$coef.live <- as.numeric(NA)
    returnVals$exp.live  <- as.numeric(NA)
    returnVals$MSE.live  <- as.numeric(NA)
    returnVals$r.live    <- as.numeric(NA)
  } else {
    coefs      <- coef(model <- nls(y ~ I(a * x^b), start = list(a = startVals, b = startVals)))
    predicted                <- coefs[1] * x^coefs[2]
    squared_error            <- (predicted - y)^2
    returnVals$coef.live  <- coefs[1]
    returnVals$exp.live   <- coefs[2]
    returnVals$MSE.live   <- mean(squared_error, na.rm = T)
    returnVals$r.live     <- 1 - (deviance(model)  / sum((y[!is.na(y)] - mean(y, na.rm = T))^2))
  }
  
  # get dead biomass parameters
  if (length(x.dead) < cutOff) {
    returnVals$coef.dead <- as.numeric(NA)
    returnVals$exp.dead  <- as.numeric(NA)
    returnVals$MSE.dead  <- as.numeric(NA)
    returnVals$r.dead    <- as.numeric(NA)
  } else {
    coefs      <- coef(model <- nls(y.dead ~ I(a * x.dead^b), start = list(a = startVals, b = startVals)))
    predicted                <- coefs[1] * x.dead^coefs[2]
    squared_error            <- (predicted - y.dead)^2
    returnVals$coef.dead  <- coefs[1]
    returnVals$exp.dead   <- coefs[2]
    returnVals$MSE.dead   <- mean(squared_error, na.rm = T)
    returnVals$r.dead     <- 1 - (deviance(model)  / sum((y.dead[!is.na(y.dead)] - mean(y.dead, na.rm = T))^2))
  }
  
  if (returnPlot == TRUE) {
    if (savePlot == TRUE) {
      fileName = paste0("Allom-", Sys.time(), ".png")
      png(filename = fileName, width = 15, height = 8, units = "cm", res = 300)
    }
    living <- c("Live", "coef.live", "exp.live", "live")
    dying <- c("Dead", "coef.dead", "exp.dead", "dead")
    
    par(mar = c(4, 4, 0.3, 0.5), fig = c(0, 0.48, 0, 1))
    biomassType <- "live"
    if (!is.na(plotTitle)) {
      title_of_plot <- paste0(plotTitle, ": live stems")
    } else {
      title_of_plot <- ""
    }
    plot(dataset[, massCol] ~ dataset[, heightCol],
         ylab = "mass (g)", xlab = "height (cm)",
         type = "n", las = 1, xlim = c(0, max(dataset[, heightCol], na.rm = T)),
         ylim = c(0, max(dataset[, massCol], na.rm = T)))
    
      points(x = dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], 
             y = dataset[, massCol][tolower(dataset[, typeCol]) %in% biomassType], 
             pch = 19, cex = 0.8, las = 1)
      # predicted values
      min.x <- min(dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], na.rm = T)
      max.x <- max(dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], na.rm = T)
      xVals <- c((min.x * 100):(max.x * 100)) / 100
      modeled <- returnVals[1, paste0("coef.", biomassType)] * xVals ^ (returnVals[1, paste0("exp.", biomassType)])
      lines(x = xVals, y = modeled, lty = i, col = "red")
      text(x = median(dataset[, heightCol], na.rm = T), y = 0.7 * max(dataset[, massCol], na.rm = T), 
           cex = 0.85, title_of_plot)
      
      biomassType <- "dead"
      if (!is.na(plotTitle)) {
        title_of_plot <- paste0(plotTitle, ": dead stems")
      } else {
        title_of_plot <- ""
      }
    par(mar = c(4, 4, 0.3, 0.5), fig = c(0.5, 1, 0, 1), new = T)
    plot(dataset[, massCol] ~ dataset[, heightCol],
         ylab = "mass (g)", xlab = "height (cm)",
         type = "n", las = 1, xlim = c(0, max(dataset[, heightCol], na.rm = T)),
         ylim = c(0, max(dataset[, massCol], na.rm = T)))
    
    points(x = dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], 
           y = dataset[, massCol][tolower(dataset[, typeCol]) %in% biomassType], 
           pch = 19, cex = 0.8, las = 1)
    # predicted values
    min.x <- min(dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], na.rm = T)
    max.x <- max(dataset[, heightCol][tolower(dataset[, typeCol]) %in% biomassType], na.rm = T)
    xVals <- c((min.x * 100):(max.x * 100)) / 100
    modeled <- returnVals[1, paste0("coef.", biomassType)] * xVals ^ (returnVals[1, paste0("exp.", biomassType)])
    lines(x = xVals, y = modeled, lty = i, col = "red")
    text(x = median(dataset[, heightCol], na.rm = T), y = 0.7 * max(dataset[, massCol], na.rm = T), 
         cex = 0.85, title_of_plot)
    if (savePlot == TRUE) {
      dev.off()
    }
  }
  
  
  returnVals
}


getAllometryParams <- function (dataset, sitesIncluded = "all", timeName = "1", siteCol = "site",
                                returnData = "FALSE", combinePlots = "FALSE", cutOff = 5, ...) {
  # function returns parameters for exponential fits (and diagnostic plots, saved to working directory)
  # for each plot. Time span used for parameterization is the entire dataset.
  
  # dataset: a dataframe with plant masses, heights etc. 
  # sitesIncluded: which sites to report data for. "all" (default), "LUM", or "TB"
  # returnData: should raw data be returned? If "TRUE", 
  #  output is a list with two elements: allometry parameters and raw data.
  # combinePlots: if `TRUE`, data from all plots are pooled (if sitesIncluded = "all", all plots are pooled). "plot" value will then include all combined plots. 
  # timeName: name of time period covered by dataset
  # cutOff: minimum number of observations to generate an allometric relationship
  
  # getParams(cwc[(cwc$season %in% "sprg 13") & (cwc$marsh %in% "LUM"), ])
  # getAllometryParams(cwc[(cwc$season %in% "sprg 13") & (cwc$marsh %in% "LUM"), ], combinePlots = "TRUE")
  
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
  
  # if necessary, clean up the dataset
  #   dataset <- mergeMonths(dataset)
  
#   # find month (robust to datasets spanning multiple months)
#   mo <- paste(paste(as.character(grep(paste(c(unique(substr(dataset$time, 1, 3))), collapse = "|"), month.abb))), 
#               collapse = " ")
#   # year reports all years covered by a dataset; specific year-month combinations are reported by moYr
#   da <- paste0(paste0("20", unique(as.character(substr(dataset$time, 5, 6)))), collapse = " ")
#   moYr <- paste0(unique(as.character(dataset$time)), collapse = " ")
#   
  # check site exclusion flag
  if ((sitesIncluded %in% "all") & (combinePlots %in% countsAsFalse)) {
    plotsInData <- unique(dataset$site)[unique(dataset$site) %in% c(TB.sites, LUM.sites)]
  } else if ((sitesIncluded %in% "all") & (combinePlots %in% countsAsTrue)) {
    # list all plots merged, if plots are being merged
    plotsInData <- paste0(unique(dataset$site)[unique(dataset$site) %in% c(TB.sites, LUM.sites)], collapse = " ")
  } else if ((sitesIncluded %in% "LUM") & (combinePlots %in% countsAsFalse)) {
    plotsInData <- unique(dataset$site)[unique(dataset$site) %in% LUM.sites]
    dataset <- dataset[dataset$site %in% LUM.sites, ]
  } else if ((sitesIncluded %in% "TB") & (combinePlots %in% countsAsFalse))  {
    plotsInData <- unique(dataset$site)[unique(dataset$site) %in% TB.sites]
    dataset <- dataset[dataset$site %in% TB.sites, ]
  } else if ((sitesIncluded %in% "LUM") & (combinePlots %in% countsAsTrue)) {
    plotsInData <- paste0(unique(dataset$site)[unique(dataset$site) %in% LUM.sites], collapse = " ")
    dataset <- dataset[dataset$site %in% LUM.sites, ]
  } else if ((sitesIncluded %in% "TB") & (combinePlots %in% countsAsTrue)) {
    plotsInData <- paste0(unique(dataset$site)[unique(dataset$site) %in% TB.sites], collapse = " ")
    dataset <- dataset[dataset$site %in% TB.sites, ]
  }
  
  if (length(plotsInData) == 0) {
    print(paste0("no plots from ", sitesIncluded, " were found. Try using `all`."))
    returnVals <- data.frame(timePeriod = timeName, 
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
    # if all sites are selected, this step begins to combine their data
    if ((combinePlots %in% countsAsTrue) & (sitesIncluded %in% "all")) {
      n <- 1
    } else if ((combinePlots %in% countsAsTrue) & (!sitesIncluded %in% "all")) {
      n <- 1
    } else {
      n <- length(plotsInData) 
    }
    
    returnVals <- data.frame(timePeriod = rep(timeName, times = n),
                             site = rep(NA, times = n), plot = rep(NA, times = n),
                             coef.live = rep(as.numeric(NA), times = n), exp.live = rep(as.numeric(NA), times = n), # y = coef * e^(exp * b)
                             MSE.live = rep(as.numeric(NA), times = n), r.live = rep(as.numeric(NA), times = n),    # diagnostics
                             coef.dead = rep(as.numeric(NA), times = n), exp.dead = rep(as.numeric(NA), times = n), 
                             MSE.dead = rep(as.numeric(NA), times = n), r.dead = rep(as.numeric(NA), times = n)
    )    
    mergeMarker <- grep("coef.live", names(returnVals))
    
    for (i in 1:n) {
      if (combinePlots %in% countsAsTrue) {
        siteName <- strsplit(plotsInData, " ")[[1]]
      }  else {
        siteName <- plotsInData[i]
      }
      
      if ((length(sitesInData) > 1) & (combinePlots %in% countsAsTrue)) {
        returnVals$site[i] <- paste0(sitesInData, collapse = " ")
      } else {
        returnVals$site[i] <- sapply(strsplit(as.character(unique(dataset$site)), "[1-9]| "), head, n = 1)[1] # old version: strsplit(as.character(unique(dataset$site)), split = c("[1-9]| "))[[i]][1] # "LUM" or "TB"
      }
      
      # `plot` reports `sitesIncluded` if data are being combined across plots
      if (combinePlots %in% countsAsFalse) {
        returnVals$plot[i]  <- plotsInData[i]
      } else if ((combinePlots %in% countsAsTrue) & (sitesIncluded %in% "all")) {
        returnVals$plot[i]  <- plotsInData[i]
      } else if ((combinePlots %in% countsAsTrue) & (!sitesIncluded %in% "all")) {
        returnVals$plot[i]  <-  returnVals$site[i]
      }
      
      ### data subset by type
      subDat       <- dataset[dataset[, siteCol] %in% siteName, ]
      returnVals[i, mergeMarker:ncol(returnVals)] <- getParams(subDat, ...)[3:10] # magic numbers. replace 3:10 with -c(1:2)
      
    }
  }
  if (returnData %in% countsAsFalse) {
    returnVals
  } else if (returnData %in% countsAsTrue) {
    returnVals <- list(parameters = returnVals, data = dataset)
  } 
}


plotAllom <- function(monthlyData, site = "LUM", type = "both", save = "TRUE", ...) {
  procData <- getAllometryParams(dataset = monthlyData, returnData = "TRUE", sitesIncluded = site, ...)
  
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
           title = paste0(procData[[1]]$timePeriod[1], " ", living[4], " stems")
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
           title = paste0(procData[[1]]$timePeriod[1], " ", dying[4], " stems")
    )
  }
  if (save %in% c("TRUE", "True", "true", "T")) {
    dev.off()
  }
}




PSC <- function(dataset, liveCol = "live", deadCol = "dead", yearCol = "year", siteCol = "site", type = "both",
                maxMinIncrement = "TRUE") {
  # peak standing crop estimates of NAPP
  # runs for entire dataset, reports results by year for each site (summary statistics)
  #   dataset = dataframe with your data
  #   liveCol = name of the column with live biomass data
  #   deadCol = name of the column with dead biomass data
  #   yearCol = name of the column with year (4 digits)
  #   siteCol = name of the column with site name
  #   type = "PSC-A" is the traditional approach, using just live material; "PSC-B" uses live + dead material in estimates of peak standing crop
  #   maxMinIncrement = indicates whether maximum-mininum biomass increment should be calculated
  
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
  countsAsTrue  <- c("T", "TRUE", "true", "True")
  
  for (h in 1:length(unique(dataset[, siteCol]))) {
    targetSite <- unique(dataset[, siteCol])[h]
    subData1 <- dataset[dataset[, siteCol] %in% targetSite, ]
    for (i in 1:length(unique(subData1[, yearCol]))) {
      targetYear <- unique(subData1[, yearCol])[i]
      subData <- subData1[subData1[, yearCol] %in% targetYear, ]
      if (type %in% "both") {
        PSC_A    <- max(subData[, liveCol], na.rm = T)
        PSC_B    <- max(subData[, liveCol], na.rm = T) + subData[, deadCol][subData[, liveCol] == PSC_A]
        if (maxMinIncrement %in% countsAsTrue) {
          maxMin_a <- PSC_A - min(subData[, liveCol], na.rm = T)
          maxMin_b <- PSC_B - min((subData[, liveCol] + subData[, deadCol]), na.rm = T)
          output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A, psc_b = PSC_B,
                               maxMin_a = maxMin_a, maxMin_b = maxMin_b)
        } else {
          output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A, psc_b = PSC_B)
        }
      } else if (type %in% "PSC-A") {
        PSC_A <- max(subData[, liveCol], na.rm = T)
        maxMin_a <- PSC_A - min(subData[, liveCol], na.rm = T)
        if (maxMinIncrement %in% countsAsTrue) {
          maxMin_a <- PSC_A - min(subData[, liveCol], na.rm = T)
          output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A, maxMin_a = maxMin_a)
        } else {
          output <- data.frame(site = targetSite, year = targetYear, psc_a = PSC_A)
        }
      } else if (type %in% "PSC-B") {
        PSC_B  <- max(subData[, liveCol], na.rm = T) + subData[, deadCol][subData[, liveCol] == max(subData[, liveCol], na.rm = T)]
        if (maxMinIncrement %in% countsAsTrue) {
          maxMin_b <- PSC_B - min((subData[, liveCol] + subData[, deadCol]), na.rm = T)
          output <- data.frame(site = targetSite, year = targetYear, psc_b = PSC_B, maxMin_b = maxMin_b)
        } else {        
          output <- data.frame(site = targetSite, year = targetYear, psc_b = PSC_B)
        }
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
                     MilnerHughes = "TRUE", EOS = "FALSE", EOS_window = 1,
                     summarize = "FALSE", timeCol = "time") {
  # requires that zoo library be loaded (time is converted to yearmon for finding EOS)
  # implements Smalley (1958) and Milner and Hughes (1968)
  # runs for entire dataset, reports results by year for each plot
  #   dataset = dataframe with your data
  #   liveCol = name of the column with live biomass data
  #   deadCol = name of the column with dead biomass data
  #   yearCol = name of the column with year (4 digits)
  #   siteCol = name of the column with plot name
  #   timeCol = column with time data (in form "%B %Y"). This only matters for the summary statistics, which report peak timing
  #   MilnerHughes      = if "TRUE", also implements Millner & Hughes 1968 (sum of positive changes in standing live biomass)
  #   summarize = if "TRUE", summary statistics (max NAPP estimates) are reported. TODO: report peak timing
  #   EOS    = "TRUE" includes a column calculating September biomass (or closest month in dataset). If there was 
  #              no sampling within some number of months (+- EOS_window) of September, value is reported as NA
  #
  # Usage examples: 
  # # Single site, single year
  #   test <- napp[(napp$site %in% "LUM1") & (napp$year %in% "2014"), 1:6]
  #   nappCalc(test)
  # # Single site, multiple years
  #   test2 <- napp[(napp$site %in% "LUM1"), 1:6]
  #   nappCalc(test2, summarize = "TRUE")[[2]]
  # # Multiple sites, multiple years
  #   test3 <- napp[, 1:6]
  #   nappCalc(test3)
  #   nappCalc(test3, summarize = "TRUE")[[2]]
  #   nappCalc(test3, summarize = "TRUE")[[2]]
  # # Last two columns are identical to 
  #   PSC(napp[, 1:6])
  
  ### error checking
  countsAsTrue  <- c("T", "TRUE", "true", "True")
  countsAsFalse <- c("F", "FALSE", "false", "False")
  
  if (!MilnerHughes %in% c(countsAsTrue, countsAsFalse)) {
    stop ("`MilnerHughes` argument isn't recognized. Input can be either `TRUE` or `FALSE`.")
  }
  if (sum(c(liveCol, deadCol, yearCol, siteCol) %in% names(dataset)) < 4) {
    stop ("Check column names. One or more column names were not found in the dataset.")
  }
  ###
  
  tempData    <- dataset
  # column names as variables:
  smalley.inc <- "smalley.inc"
  smalley     <- "smalley"
  live.inc    <- "live.inc" # VTS1975's delL
  dead.inc    <- "dead.inc" # VTS1975's delD
  MH          <- "MH"
  maxMin      <- "maxMin"
  # Valiela, Teal, Sass 1975
  eV          <- "VTS1975.inc"
  VTS         <- "VTS1975"
  PSC_A       <- "psc.live"
  PSC_B       <- "psc.tot"
  EOS_col     <- "eos"
  
  # define acceptable window for EOS measurement
  EOS_target <- as.numeric(as.yearmon("Sep 2016", format = "%B")) - 2016
  EOS_high   <- as.numeric(as.yearmon(paste0(month.abb[grep("Sep", month.abb) + EOS_window], "2016"), format = "%B")) - 2016
  EOS_low   <- as.numeric(as.yearmon(paste0(month.abb[grep("Sep", month.abb) - EOS_window], "2016"), format = "%B")) - 2016
  
  # more variables than necessary are appended to dataset
  tempData[, EOS_col] <- tempData[, PSC_B] <- tempData[, PSC_A] <- tempData[, VTS] <- tempData[, MH] <- tempData[, smalley] <- tempData[, smalley.inc] <- 
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
      
      # sum positive *live* biomass increments for Millner & Hughes 1968
      subData2[, "MH.inc"] <- subData2[, live.inc]
      # use only positive increments to calc MH NAPP
      subData2[, "MH.inc"][subData2[, "MH.inc"] < 0] <- 0
      subData2[, MH][!is.na(subData2[, "MH.inc"])] <- cumsum(subData2[, "MH.inc"][!is.na(subData2[, "MH.inc"])])
      
      # sum e from Valiela, Teal, Sass 1975
      subData2[, VTS][!is.na(subData2[, eV])] <- cumsum(subData2[, eV][!is.na(subData2[, eV])])
      
      # sum peak standing crop increments
      subData2[, PSC_A][!is.na(subData2[, liveCol])]                       <- cummax(subData2[, liveCol][!is.na(subData2[, liveCol])])
      # PSC_B reflects maximum summed biomass when live biomass is at its peak. 
      subData2[, PSC_B][!is.na(subData2[, liveCol] + subData2[, deadCol])] <- max(subData2[, liveCol][!is.na(subData2[, liveCol])]) + subData2[, deadCol][subData2[, liveCol][!is.na(subData2[, liveCol])] == max(subData2[, liveCol][!is.na(subData2[, liveCol])])]
      
      # find "end of season" biomass
      # may need to ensure that monthly data exist? could pose a problem when assigning EOS_biomass and querying months.
      if(EOS %in% countsAsTrue) {
        # if no sampling occurred in September, widen window
        if (round(EOS_target, 2) %in% round((as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol])))), 2)) {
          EOS_biomass <- subData2[, liveCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) == EOS_target)] +
            subData2[, deadCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) == EOS_target)]
        # } else if(sum(c(round(seq(from = EOS_low, to = EOS_high, by = 1/12), 2) %in% round((as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol])))), 2))) >= 1) {
          # EOS_biomass <- max(subData2[, liveCol][(as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) >= EOS_low) | (as.numeric(as.yearmon(subData2[, timeCol])) - floor(as.numeric(as.yearmon(subData2[, timeCol]))) <= EOS_high)], na.rm = TRUE)
        } else {
          EOS_biomass <- NA
        }
        subData2[, EOS_col][subData2[, liveCol] + subData2[, deadCol] == EOS_biomass] <- EOS_biomass
      }
      
      
      # add to output dataframe
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), smalley]  <- subData2[, smalley]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), MH]       <- subData2[, MH]
      # MH should resolve to max-min if first and final biomass data = 0
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), maxMin]   <- max(subData2[, liveCol], na.rm = T) - min(subData2[, liveCol], na.rm = T)
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), VTS]      <- subData2[, VTS]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_A]    <- subData2[, PSC_A]
      tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), PSC_B]    <- subData2[, PSC_B]
      if(EOS %in% countsAsTrue) {
        tempData[(tempData[, siteCol] %in% targetSite) & (tempData[, yearCol] %in% targetYear), EOS_col]   <- subData2[, EOS_col]
      }
    }
  }
  
  if (summarize %in% countsAsFalse) {
    output <- tempData
  } else if (summarize %in% countsAsTrue) {
    for (i in 1:length(unique(tempData[, siteCol]))) {
      targetSite <- unique(tempData[, siteCol])[i]
      subData1 <- tempData[tempData[, siteCol] %in% targetSite, ]
      for (j in 1:length(unique(subData1[, yearCol]))) {
        targetYear <- unique(subData1[, yearCol])[j]
        subData2   <- subData1[subData1[, yearCol] %in% targetYear, ] 
        
        intData <- data.frame(
          site = targetSite, 
          year = targetYear,
          mean.live    = mean(subData2[, liveCol], na.rm = T),
          mean.dead    = mean(subData2[, deadCol], na.rm = T),
          napp.smalley = max(subData2[, smalley], na.rm = T), 
          napp.MH      = max(subData2[, MH], na.rm = T), 
          napp.maxMin  = max(subData2[, maxMin], na.rm = T), 
          napp.VTS     = max(subData2[, VTS], na.rm = T), 
          napp.psc.a   = max(subData2[, PSC_A], na.rm = T), 
          napp.psc.b   = max(subData2[, PSC_B], na.rm = T), 
          n            = sum(!is.na(subData2[, liveCol])),
          t.smalley    = ifelse(is.finite(max(subData2[, smalley], na.rm = T)), as.character(subData2[, timeCol][which.max(subData2[, smalley])]), NA),
          t.MH         = ifelse(is.finite(max(subData2[, MH], na.rm = T)), as.character(subData2[, timeCol][which.max(subData2[, MH])]), NA),
          t.vts        = as.character(subData2[, timeCol][which.max(subData2[, VTS])]),
          t.psc.a      = as.character(subData2[, timeCol][which.max(subData2[, PSC_A])]),
          t.psc.b      = as.character(subData2[, timeCol][which.max(subData2[, PSC_B])])
        )
        
        # ADD EOS to intdata if EOS == TRUE
        if (EOS %in% countsAsTrue) {
          intData$napp.EOS <- ifelse(sum(is.na(subData2[, EOS_col])) == length(subData2[, EOS_col]), 
                                      NA, max(subData2[, EOS_col], na.rm = T)) 
        }
        
        
        if (i != 1 ) {
          finalData <- rbind(finalData, intData)
        } else if ((i == 1) & (j != 1)) {
          finalData <- rbind(finalData, intData)
        } else {
          finalData <- intData
        }
        
      }
    }
    finalData$napp.MH[!is.finite(finalData$napp.MH)] <- NA
    output <- list(intervalData = tempData, summary = finalData)
  }
  output
}


# note 150907: changed to use y = a*x^b
predictBiomass <- function(plotData = cwc, monthYear, plot, quadrat = 0.25, 
                           removeTrainingData = "FALSE", returnData = "TRUE", start_nls = 0.03,
                           coefReturn = "FALSE", cutoff = 5, moYrFormat = "%B-%y", 
                           typeCol = "type", timeCol = "monthYear", heightCol = "hgt",
                           siteCol = "site" , massCol = "mass") {
  # function to apply allometry to other time points (same plot). Generates a residual (observed - predicted)
  # allometryData:  object with allometry parameters (cwc.params)
  # plotData:       object with plot data on plant heights and masses (cwc)
  # summaryData:    data with total biomass per plot, per time interval (maybe not needed)
  # monthYear:      month and year to use for estimating biomass. monthYear should be of the form, e.g., "May-13" for May 2013
  # plot:           which plot's data will be used
  # removeTrainingData: determines whether training data should also be modeled
  # returnData:     indicates whether predicted is returned by the function
  #   
  #   plotData      <- cwc
  #   summaryData   <- CWC.plots # plot.totals has combined live/dead, CWC.plots distinguishes live/dead
  
  # TODO: add site name to coef dataframe
  
  countsAsTrue <- c("TRUE", "true", "T", "True")
  plotSize     <- quadrat^2
  
  # isolate plot x's data (and exclude month(s) used to parameterize model, if desired)
  if (removeTrainingData %in% countsAsTrue) {
    trainingDataSubset <- plotData[(as.character(plotData[, siteCol]) %in% plot) & (plotData[, timeCol] %in% monthYear), ]
    # newData: full dataset, modified before being output. 'removeTrainingData' removes the training data from this dataset.
    newData            <- plotData[(!as.character(plotData[, siteCol]) %in% plot) & (!plotData[, timeCol] %in% monthYear), ]
  } else {
    trainingDataSubset <- plotData[(as.character(plotData[, siteCol]) %in% plot) & (plotData[, timeCol] %in% monthYear), ]
    newData            <- plotData
  }
  
  test.live <- nrow(trainingDataSubset[(trainingDataSubset[, typeCol] %in% "LIVE") & (trainingDataSubset[, timeCol]  %in% monthYear), ]) > cutoff
  test.dead <- nrow(trainingDataSubset[(trainingDataSubset[, typeCol] %in% "DEAD") & (trainingDataSubset[, timeCol]  %in% monthYear), ]) > cutoff
  
  # generate allometry params (allows multiple months to be combined)
  if (test.live) {
    y.live <- trainingDataSubset[(trainingDataSubset[, timeCol] %in% monthYear) & (trainingDataSubset[, typeCol] %in% "LIVE"), massCol]
    x.live <- trainingDataSubset[(trainingDataSubset[, timeCol] %in% monthYear) & (trainingDataSubset[, typeCol] %in% "LIVE"), heightCol]
    y.live2 <- y.live[!is.na(y.live) & !is.na(x.live)]
    x.live2 <- x.live[!is.na(y.live) & !is.na(x.live)]
    live.coefs <- coef(model <- nls(y.live2 ~ I(a * x.live2^b), start = list(a = start_nls, 
                                                                             b = start_nls)))
  } else {
    live.coefs <- c(NA, NA)
  }
  # force through zero?
  #   plot(y.live2 ~ x.live2)
  #   lines(1:150, y = live.coefs[1] * exp(live.coefs[2] * 1:150))
  
  if (test.dead) {
    y.dead <- trainingDataSubset[(trainingDataSubset[, timeCol]  %in% monthYear) & (trainingDataSubset[, typeCol] %in% "DEAD"), massCol]
    x.dead <- trainingDataSubset[(trainingDataSubset[, timeCol]  %in% monthYear) & (trainingDataSubset[, typeCol] %in% "DEAD"), heightCol]
    y.dead2 <- y.dead[!is.na(y.dead) & !is.na(x.dead)]
    x.dead2 <- x.dead[!is.na(y.dead) & !is.na(x.dead)]
    dead.coefs <- coef(model <- nls(y.dead2 ~ I(a * x.dead2^b), start = list(a = start_nls, 
                                                                             b = start_nls)))
  } else {
    dead.coefs <- c(NA, NA)
  }
  #   plot(y.dead2 ~ x.dead2)
  #   lines(1:150, y = dead.coefs[1] * exp(dead.coefs[2] * 1:150))
  
  a.live <- live.coefs[1]
  b.live <- live.coefs[2]
  a.dead <- dead.coefs[1]    # coef
  b.dead <- dead.coefs[2]
  
  # apply allometry from t1 to other time points
  newData$predicted <- NA
  newData$predicted[newData[, typeCol] %in% "LIVE"] <- a.live * newData[newData[, typeCol] %in% "LIVE", heightCol]^b.live
  newData$predicted[newData[, typeCol] %in% "DEAD"] <- a.dead * newData[newData[, typeCol] %in% "DEAD", heightCol]^b.dead
  # error per stem (difference as percent of observed)
  newData$stem.err <- (newData$predicted - newData[, massCol]) / newData[, massCol]
  
  
  ##########
  # for each site-month-type combination, calculate error metrics
  # after checking to see if observations exceed cutoff 
  for (i in 1:length(unique(newData[, timeCol]))) { # for each month
    targetTime <- unique(newData[, timeCol])[i]
    subData    <- newData[newData[, timeCol] %in% targetTime, ]
    
    for (j in 1:length(unique(subData[, siteCol]))) { # for each site
      targetSite <- unique(subData[, siteCol])[j]
      subData2   <- subData[subData[, siteCol] %in% targetSite, ]
      
      for (k in 1:length(unique(subData2[, typeCol]))) { # for live/dead
        targetType <- unique(subData2[, typeCol])[k]
        subData3   <- subData2[subData2[, typeCol] %in% targetType, ]
        
        ### perform main routine ###
        if(!is.null(subData3)) {
          intOutput_LD <- data.frame(
            plot      = targetSite,
            monthYear = targetTime,
            type      = targetType,
            stem.err  = median(subData3$stem.err, na.rm = T),
            # error per plot
            mass.obs  = sum(subData3[, massCol], na.rm = T)/ plotSize,
            mass.pred = sum(subData3$predicted, na.rm = T) / plotSize,
            err.pct   = (sum(subData3[, massCol], na.rm = T) - sum(subData3$predicted, na.rm = T) ) / sum(subData3[, massCol], na.rm = T)
          )
        } else if (is.null(subData3)) {
          intOutput_LD <- data.frame(
            plot      = targetSite,
            monthYear = targetTime,
            type      = targetType,
            stem.err  = as.numeric(NA),
            mass.obs  = as.numeric(NA),
            mass.pred = as.numeric(NA),
            err.pct   = as.numeric(NA)
          )
        }
        
        ### now, merge ###
        if (k == 1) {
          output_LD <- intOutput_LD
        } else {
          output_LD <- rbind(output_LD, intOutput_LD)
        }
      }
      if (j == 1) {
        output_site <- output_LD
      } else {
        output_site <- rbind(output_site, output_LD)
      }
    }
    if (i == 1) {
      output <- output_site
    } else {
      output <- rbind(output, output_site)
    }
  }
  
  ##########
  ### Process output a little bit
  ### make sure that an absence of biomass has not been interpreted as an absence of sampling
#   nrow(output) # up to 207
  for (i in 1:length(unique(output$plot))) {
    for (j in 1:length(unique(output$monthYear))) {
      plot    <- unique(output$plot)[i]
      time    <- unique(output$monthYear)[j]
      subData <- output[(output$plot %in% plot) & (output$monthYear %in% time), ]
      if (nrow(subData) == 1) {
        fillData                     <- subData[1, ]
        fillData$type                <- ifelse(fillData$type %in% "LIVE", "DEAD", "LIVE")
        fillData[, 4:ncol(fillData)] <- as.numeric(0)
        output <- rbind(output, fillData)
      }
    }
    rownames(output) <- 1:nrow(output)
  }
#   nrow(output) # up to 5 larger (212)
  
  output$year <- paste0("20", substr(as.character(output$monthYear), 5, 6))
  output <- output[order(output$year, as.yearmon(output$monthYear, moYrFormat)), ]
  ##########
  
  
  if (returnData %in% countsAsTrue) {
    if (coefReturn %in% countsAsTrue) {
      coef.df <- data.frame(coef.live = a.live,
                            exp.live  = b.live,
                            coef.dead = a.dead,
                            exp.dead  = b.dead
      )
      output <- list(summary = output, data = newData, coefs = coef.df)
    } else {
      output <- list(summary = output, data = newData)
    }
  } else {
    if (coefReturn %in% countsAsTrue) {
      coef.df <- data.frame(coef.live = a.live,
                            exp.live  = b.live,
                            coef.dead = a.dead,
                            exp.dead  = b.dead
      )
      output <- list(summary = output, coefs = coef.df)
    } else {
      output
    }
  }
}




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

seasonLabel <- function(data, monthYearColumn = "moYr", siteColumn = "site", seasons = NA, year = c(13:15)) {
  # labels each row's season. 'monthYearColumn' should be of the form %b-%y
  
  if(is.na(seasons)) {
    seasons <- list(
      # spring: Mar Apr May
      sprg = c("Mar", "Apr", "May"),
      # summer: Jun Jul Aug
      sumr = c("Jun", "Jul", "Aug"),
      # fall: Sep Oct Nov
      fall = c("Sep", "Oct", "Nov"),
      # winter: Dec Jan Feb
      wint = c("Dec", "Jan", "Feb")
    )
  }
  
  data[, "season"] <- as.character(NA)
  
  for (h in 1:length(unique(data[, siteColumn]))) {
    targSite <- as.character(unique(data[, siteColumn])[h])
    for (i in 1:length(seasons)) {
      for (j in 1:length(year)) {
        # account for winter spanning two years
        if (i == length(seasons)) {
          targetDates <- paste0(seasons[[i]], "-", c(year[j], year[j] + 1, year[j] + 1))
        } else {
          targetDates <- paste0(seasons[[i]], "-", year[j])
        }  
        data[, "season"][(data[, monthYearColumn] %in% targetDates) & (data[, siteColumn] %in% targSite)] <- paste(names(seasons)[i], year[j])
      }
    }
  }
  invisible(data)
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
    "smalley" = "Smalley 1959", 
    "MH" = "Milner Hughes 1968", 
    "maxMin" = "Max-Min",
    "VTS" = "Valiela et al. 1975", 
    "psc.live" = "Peak (live)",
    "psc.tot" = "Peak (live + dead)",
    "eos"     = "End-of-season live"
  )
  return(names_li[value])
}


allomLabeller1 <- function(variable, value){
  value <- droplevels(value)
  names_li <- list(
    "live.mean" = "Live stems", 
    "dead.mean" = "Dead stems",
    
    "sumr" = "summer",
    "wint" = "winter",
    "fall" = "Fall", 
    "sprg" = "Spring"#,
#     "LUM" = "LUMCON",
#     "TB-B" = "Lake Barre",
#     "TB-A" = "Bay La Fleur"
  )
  return(names_li[value])
}


# 
# chlLabeller <- function(variable, value){
#   value <- droplevels(value)
#   names_li <- list(
#     "chlA_ugcm2" = expression("Chl. a ("*mu*"g"%.%"cm"^2*")"), 
#     "pgmt_ugcm2" = expression("pigment ("*mu*"g"%.%"cm"^2*")"),
#     "chlA_ugg"   = expression("Chl. a ("*mu*"g"%.%"g"^-1*")"),
#     "pgmt_ugg"   = expression("pigment ("*mu*"g"%.%"g"^-1*")")
#   )
#   return(names_li[value])
# }

chlLabeller <- function(variable, value){
  value <- droplevels(value)
  names_li <- list(
    "chlA_inv" = expression("Chl. a (mg"%.%"m"^2*")"), 
    "pgmt_inv" = expression("pigment (mg"%.%"m"^2*")"),
    "chlA_conc"   = expression("Chl. a (mg"%.%"g"^-1*")"),
    "pgmt_conc"   = expression("pigment (mg"%.%"g"^-1*")")
  )
  return(names_li[value])
}

zeroToNA <- function(dataset, cols = c(1:ncol(dataset))) {
  # function converts all zeroes in a dataset to NAs
  for (i in 1:length(cols)) {
    dataset[, cols[i]][dataset[, cols[i]] == 0] <- NA
  }
  dataset
}



### a function to make correlation matrices with significance stars
### from http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html
### usage: kable(corstarsl(as.matrix(data))) 
corstarsl <- function(x){
  require(Hmisc)
  x <- as.matrix(x)
  R <- rcorr(x)$r
  p <- rcorr(x)$P
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
} 


### A function to create a unique ID column in a dataset, based on combinations of input columns
createUniqueID <- function (dataset, inputColNames) {
  n <- length(inputColNames)
  ID_col <- as.character(dataset[, inputColNames[1]])
  for (i in 2:n) {
    ID_col <- paste0(ID_col, "-", as.character(dataset[, inputColNames[i]]))
  }
  ID_col
}





### function creates weights for multcomp::nls
### from: https://rmazing.wordpress.com/2012/07/19/a-weighting-function-for-nls-nlslm/
wfct <- function(expr)
{
  expr <- deparse(substitute(expr))
  
  ## create new environment
  newEnv <- new.env()
  
  ## get call
  mc <- sys.calls()[[1]]
  mcL <- as.list(mc)
  
  ## get data and write to newEnv
  DATA <- mcL[["data"]]
  DATA <- eval(DATA)
  DATA <- as.list(DATA)
  NAMES <- names(DATA)
  for (i in 1:length(DATA)) assign(NAMES[i], DATA[[i]], envir = newEnv)
  
  ## get parameter, response and predictor names
  formula <- as.formula(mcL[[2]])
  VARS <- all.vars(formula)
  RESP <- VARS[1]
  RHS <- VARS[-1]
  PRED <- match(RHS, names(DATA))
  PRED <- names(DATA)[na.omit(PRED)]
  
  ## calculate variances for response values if "error" is in expression
  ## and write to newEnv
  if (length(grep("error", expr)) > 0) {
    y <- DATA[[RESP]]
    x <- DATA[[PRED]]
    ## test for replication
    if (!any(duplicated(x))) stop("No replicates available to calculate error from!")
    ## calculate error
    error <- tapply(y, x, function(e) var(e, na.rm = TRUE))
    error <- as.numeric(sqrt(error))
    ## convert to original repititions
    error <- rep(error, as.numeric(table(x)))
    assign("error", error, envir = newEnv)
  }
  
  ## calculate fitted or residual values if "fitted"/"resid" is in expression
  ## and write to newEnv
  if (length(grep("fitted", expr)) > 0 || length(grep("resid", expr)) > 0) {
    mc2 <- mc
    mc2$weights <- NULL
    MODEL <- eval(mc2)
    fitted <- fitted(MODEL)
    resid <- residuals(MODEL)
    assign("fitted", fitted, newEnv)
    assign("resid", resid, newEnv)
  }
  
  ## return evaluation in newEnv: vector of weights
  OUT <- eval(parse(text = expr), envir = newEnv)
  return(OUT)
}
