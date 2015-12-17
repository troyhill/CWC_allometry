### script re-calculates NAPP data from literature


nappLit <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_literatureNAPP_151105.txt", skip = 1)

# remove all of the Hopkinson et al. 1978 references (use 1980).
nappLit <- nappLit[!nappLit$source %in% "Hopkinson et al. 1978", ]
rownames(nappLit) <- 1:nrow(nappLit)

kable(corstarsl(as.matrix(nappLit[c(1:5, 6, 9:15, 17:20), c(12, 11, 10, 15, 8, 9)]))) # includes Gulf coast, Georgia and North Carolina
kable(corstarsl(as.matrix(nappLit[c(1:4, 6, 9:14, 17:20), c(12, 11, 10, 15, 8, 9)]))) # without north carolina, nothing significant
plot(nappLit[c(1:4, 6, 9, 17:20), c(12, 11, 10, 15, 8, 9)])

kable(corstarsl(as.matrix(nappLit[c(1:4, 6, 9, 17:20), c(12, 11, 10, 15, 8, 9)]))) # just LA, MS

salit <- read.delim("C:/RDATA/SPAL_allometry/data_LUM123/data_literatureValues_150930_.txt")
salit$moYr <- as.yearmon(paste0(substr(salit$month, 1, 3), "-", salit$year), "%b-%Y")
salit$julian <- julian(as.POSIXct(salit$moYr), origin = as.POSIXct("1970-01-01", tz = "GMT"))

# calculate julian days for White et al. 1978 (sampled at 3-week intervals)
startdate <- as.Date(paste0(salit$year[!is.na(salit$day)][1], "-", substr(salit$month[!is.na(salit$day)][1], 1, 3),"-", salit$day[!is.na(salit$day)][1]), "%Y-%b-%d")
for (i in 1:sum(!is.na(salit$day))) {
  salit$julian[!is.na(salit$day)][i] <- julian(as.POSIXct(startdate), origin = as.POSIXct("1970-01-01", tz = "GMT")) + salit$day[!is.na(salit$day)][i] - salit$day[!is.na(salit$day)][1]
}

shortOnly <- salit[-c(grep("streamside", salit$location)), ]

plot(salit$live[salit$state %in% c("LA", "MS", "AL")] ~ 
       salit$julian[salit$state %in% c("LA", "MS", "AL")], 
       pch = 19, cex = 0.7, las = 1, col = "red",
       ylab = "monthly live biomass along Gulf Coast (g/m2)", xlab = "days since 1970")
points(shortOnly$live[shortOnly$state %in% c("LA", "MS", "AL")] ~ 
         shortOnly$julian[shortOnly$state %in% c("LA", "MS", "AL")], 
       pch = 19, cex = 0.7)

y.err.pos <- ifelse(!is.na(salit$live.se[salit$state %in% c("LA", "MS", "AL")]),  salit$live[salit$state %in% c("LA", "MS", "AL")] + salit$live.se[salit$state %in% c("LA", "MS", "AL")], salit$live[salit$state %in% c("LA", "MS", "AL")])
y.err.neg <- ifelse(!is.na(salit$live.se[salit$state %in% c("LA", "MS", "AL")]),  salit$live[salit$state %in% c("LA", "MS", "AL")] - salit$live.se[salit$state %in% c("LA", "MS", "AL")], salit$live[salit$state %in% c("LA", "MS", "AL")])

segments(salit$julian[salit$state %in% c("LA", "MS", "AL")],
         y.err.pos,
         salit$julian[salit$state %in% c("LA", "MS", "AL")],
         y.err.neg
         )

# how to convert julian days to date class?


### calculate productivity using literature
c(grep("Darby", salit$source) & which(salit$year %in% c(2004, 2005)))

grep("Darby", salit$source)
max(salit[which(salit$year %in% c(2004, 2005)), "live"]) - min(salit[which(salit$year %in% c(2004, 2005)), "live"])
nappCalc(salit[which(salit$year %in% c(2004, 2005)), ], siteCol = "location", timeCol = "moYr", summarize = "TRUE")

napVals <- nappCalc(salit[!is.na(salit$dead), ], siteCol = "location", timeCol = "moYr", summarize = "TRUE")$summary

names(napVals)
plot(napVals$napp.VTS ~ napVals$napp.smalley)
summary(lm(napVals$napp.VTS ~ napVals$napp.smalley))
