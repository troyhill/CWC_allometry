library(plyr)


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
  sitesInData <- levels(as.factor(siteOnly)) 

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

    #     plotDate   <- paste0(siteName, i, "-", moYr)
    x          <- dataset$hgt[(dataset$type %in% c("Live", "LIVE", "live")) & (dataset$site %in% siteName)]
    y          <- dataset$mass[(dataset$type %in% c("Live", "LIVE", "live")) & (dataset$site %in% siteName)]
    coefs      <- coef(model <- nls(y ~ I(a * exp(b * x)), start = list(a = 0.1, b = 0.1)))
    predicted           <- coefs[1] * exp(coefs[2] * x)
    squared_error       <- (predicted - y)^2
    
    # save relevant data
#     if (combinePlots %in% countsAsTrue) {
#       returnVals$plot[i] <- plotsInData
#     }  else {
#       returnVals$plot[i] <- siteName
#     }
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

    returnVals$coef.live[i] <- coefs[1]
    returnVals$exp.live[i]  <- coefs[2]
    returnVals$MSE.live[i]  <- sum(squared_error, na.rm = T)
    returnVals$r.live[i]    <- sqrt(1 - (deviance(model)  / sum((y[!is.na(y)] - mean(y, na.rm = T))^2)))
    
    
    # do the same for dead stems, if there's more than two stems
    if (length(dataset$hgt[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]) > 2) {
        x.dead     <- dataset$hgt[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]
        y.dead     <- dataset$mass[(dataset$type %in% c("Dead", "DEAD", "dead")) & (dataset$site %in% siteName)]
        coefs.dead <- coef(model.dead <- nls(y.dead ~ I(a * exp(b * x.dead)), start = list(a = 0.1, b = 0.1)))
        
        predicted.dead      <- coefs.dead[1] * exp(coefs.dead[2] * x.dead)
        squared_error.dead  <- (predicted.dead - y.dead)^2
        
        
        returnVals$coef.dead[i] <- coefs.dead[1]
        returnVals$exp.dead[i]  <- coefs.dead[2]
        returnVals$MSE.dead[i]  <- sum(squared_error.dead, na.rm = T)
        returnVals$r.dead[i]    <- sqrt(1 - (deviance(model.dead)  / sum((y.dead[!is.na(y.dead)] - mean(y.dead, na.rm = T))^2)))
    } 
    
#     Leaving this code here until we decide whether to use R or Excel models
# 
#     oldParams <- read.delim(ExcelParameterFile, skip = 1)
#     colNos <- list(c(2, 3), c(4, 5), c(6, 7))
#     # find Excel coefficients by date matching
#     returnVals[paste0(siteName, i, ".Excel.coef")] <- Excel.coef <- oldParams[oldParams$Month %in% moYr, colNos[[i]][1]]
#     returnVals[paste0(siteName, i, ".Excel.exp")]  <- Excel.exp  <-oldParams[oldParams$Month %in% moYr, colNos[[i]][2]]
#     
#     # calculate Excel diagnostics
#     E.predicted     <- Excel.coef * exp(Excel.exp * x)
#     E.squared_error <- (E.predicted - y)^2
#     returnVals[paste0(siteName, i, ".Excel.MSE")]  <- sum(E.squared_error,  na.rm = T)
#     returnVals[paste0(siteName, i, ".Excel.corr")] <- cor(y[!is.na(y)], E.predicted[!is.na(y)])
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


