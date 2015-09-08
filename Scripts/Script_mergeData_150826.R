###############################
### Processing LUMCON CWC sites
###############################

# This script processes the aboveground plant data for the SPAL allometry project. 
# Author: Troy Hill, thill@lumcon.edu

# Notes:
# the 2013 data are organized by plot, whereas the 2014/2015 data are organized by month.
# I exclude "bag scraps" and downed litter from this analysis. Perhaps it can be incorporated in the 
# future

# install.packages(c("zoo", "plyr")) # if 'zoo' and 'plyr' aren't installed.
library(zoo)
library(plyr)

# load custom functions
source("C:/RDATA/SPAL_allometry/CWC_allometry/CWC_functions.R") # see https://github.com/troyhill/CWC_allometry/blob/master/CWC_functions.R for updated version



#####
### process data collected since March 2014 (organized by month)
#####

# load data
directory <- "C:/RDATA/SPAL_allometry/data_LUM123/"
filenames <- paste0(directory, list.files(directory, pattern = "data_CWC_20"))
dat <- lapply(filenames, read.delim)

# Calculate allometric parameters and data
data14 <- batch(inputList = dat, fun = getAllometryParams, returnData = "TRUE")


#####
### process data collected prior to March 2014 (organized by site)
#####

# load data 
directory <- "C:/RDATA/SPAL_allometry/data_LUM123/"
filenames <- paste0(directory, list.files(directory, pattern = "2013.txt"))
rawData2013 <- lapply(filenames, read.delim, skip = 3) 

# get plot names (this isn't in all data files)
sites <- substr(filenames, 42, nchar(filenames) - 9)


# combine files into a single dataframe
for(i in 1:length(rawData2013)) {
  # add site name
  inputData        <- rawData2013[[i]][, c(1:2, 4:5, 7:9, 11)]
  names(inputData) <- c("site", "time", "ID", "type", "hgt", "tin", "tin_plant", "mass")
  inputData[, "site"] <- sites[i]
  if (i != 1) {
    data13 <- rbind(data13,  inputData)
  } else {
    data13 <-  inputData
  }
}

# remove blank lines from data
data13 <- data13[!data13$ID %in% "", ]

### re-label dates to match 2014/15 data
dateList <- strsplit(x = as.character(data13$time), split = "/")
for (i in 1:length(dateList)) {
  month <- month.abb[as.numeric(dateList[[i]][1])]
 if (i != 1) {
   dateList.proc <- mapply(rbind, dateList.proc, dateList[[i]][c(1, 3)], SIMPLIFY = FALSE)
 } else {
   dateList.proc <- dateList[[i]][c(1, 3)]
 }
}

data13$time <- as.character(paste0(month.abb[as.numeric(dateList.proc[[1]])], "-", substr(as.character(dateList.proc[[2]]), 3, 4)))
data13$ID   <- as.character(data13$ID)



#####
### Merge raw data, add/change some columns
#####
cwc <- rbind(data13, data14[[2]])

# add marsh name
cwc <- marshName(cwc)

# samples to check:
cwc[(cwc$monthYear %in% "Jun-14") & (cwc$hgt < 5), ] # very odd; tiny stem, large mass tin 2131
cwc[(cwc$monthYear %in% "Jun-14") & (cwc$hgt < 5), c("mass", "hgt")] <- NA # remove for now
cwc[cwc$monthYear %in% "Jan-14",]                    # very noisy data
plot(cwc$mass[cwc$monthYear %in% "Jun-14"] ~ cwc$hgt[cwc$monthYear %in% "Jun-14"])
plot(cwc$mass[cwc$monthYear %in% "Jan-14"] ~ cwc$hgt[cwc$monthYear %in% "Jan-14"])
#####


### remove "bag scrap" samples
cwc           <- cwc[!cwc$ID %in% "bag scrap", ]
cwc$type[cwc$type %in% c("Live")] <- "LIVE"
cwc$type[cwc$type %in% c("Dead")] <- "DEAD"
cwc           <- droplevels(cwc)
cwc$monthYear <- cwc$time
cwc$time      <- as.yearmon(cwc$time, "%b-%y")
cwc$type      <- as.character(cwc$type)



#####
### get allometry for early period, merge allometry datafiles
#####
months <- unique(data13$time)
for (i in 1:length(months)) {
  sub <- data13[data13$time %in% months[i], c("site", "time", "type", "ID", "hgt", "tin", "tin_plant", "mass")]
  params2013_temp <- getAllometryParams(dataset = sub)
  if (i != 1) {
    params2013 <- rbind(params2013, params2013_temp)
  } else {
    params2013 <- params2013_temp
  }
}

# merged allometric parameters
cwc.params <- rbind(params2013, data14[[1]]) 





### save the data as .csv files
# write.csv(cwc.params, file = "C:/RDATA/SPAL_allometry/CWC_parameters_150825.csv")
# write.csv(cwc, file = "C:/RDATA/SPAL_allometry/CWC_allometryData_150825.csv")
