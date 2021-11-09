##------------------------------------------------------------------------------##
## File:        prepare.data.R
##
## Description: This file (1) compiles and (2) prepares the data used in LW.
##
##------------------------------------------------------------------------------##
rm(list = ls())
source("utilities.R")

## Load time series library
if (!require("tis")) {install.packages("tis"); library('tis')}
## Load library to read excel files
if (!require("xlsx")) {install.packages("xlsx"); library('xlsx')}

data.out.start <- c(1948,1)

##------------------------------------------------------------------------------##
## Get Raw Data
##------------------------------------------------------------------------------##

rstar.data <- read.xlsx("rstar_data.xlsx", sheetName="Quarterly", rowIndex=2:4000, colIndex=1:10, as.data.frame=TRUE, header=FALSE)
names(rstar.data) <- c("date","gdp","pce.index","core.pce.index","imported.oil","core.imports","frbus.oil","frbus.imports","frbny.ffr","ffr")

## Set the start and end dates of the data used in the estimation
data.start <- c(year(rstar.data$date[1]), quarter(rstar.data$date[1]))
data.end   <- c(year(rstar.data$date[length(rstar.data$date)]), quarter(rstar.data$date[length(rstar.data$date)]))

## Define variables as tis objects
for (var in names(rstar.data)[2:length(names(rstar.data))]) {
    assign(var, tis(rstar.data[,var], start = data.start, tif = 'quarterly'))
}


##------------------------------------------------------------------------------##
## Prepare Data
##------------------------------------------------------------------------------##

## Take log of real GDP
gdp.log <- log(gdp)

## Splice frbus.imports with core.imports at 1967:Q1
## Beginning of data - 1966:Q4 : PMO (frbus.imports) * (core.imports[1967:Q1] / frbus.imports[1967:Q1])
## 1967:Q1 - End of data : core.imports
scaled.frbus.imports <- window(frbus.imports, end = c(1966,4))
scaled.frbus.imports <- scaled.frbus.imports * (core.imports[ti(c(1967,1),'quarterly')] / frbus.imports[ti(c(1967,1),'quarterly')])
import.index <- mergeSeries(scaled.frbus.imports,
                            window(core.imports, start=c(1967,1)))

## Import price inflation (annualized)
import.price.inflation <- 400*log(import.index/Lag(import.index, k=1))

## Splice frbus.oil with imported.oil at 1967:Q1
## Beginning of data - 1966:Q4 : POIL (frbus.oil) * (imported.oil[1967:Q1] / frbus.oil[1967:Q1])
## 1967:Q1 - End of data : imported.oil
scaled.frbus.oil <- window(frbus.oil, end = c(1966,4))
scaled.frbus.oil <- scaled.frbus.oil * (imported.oil[ti(c(1967,1),'quarterly')] / frbus.oil[ti(c(1967,1),'quarterly')])
oil.index <- mergeSeries(scaled.frbus.oil,
                         window(imported.oil, start=c(1967,1)))

## Oil price inflation (annualized)
oil.price.inflation <- 400*log(oil.index/Lag(oil.index, k=1))

## Create annualized inflation series using the price indexes
pce.inflation      <- 400*log(pce.index/Lag(pce.index, k=1))
core.pce.inflation <- 400*log(core.pce.index/Lag(core.pce.index, k=1))

## Splice pce.inflation and core.pce.inflation in 1959:Q2
inflation <- mergeSeries(window(pce.inflation, end = c(1959,1)),
                         window(core.pce.inflation, start = c(1959,2)))

## Compute inflation expectations series
y.pi <- (inflation + Lag(inflation,k=1) + Lag(inflation,k=2) + Lag(inflation,k=3))/4
x.pi <- cbind(Lag(inflation,k=4),
              Lag(inflation,k=5),
              Lag(inflation,k=6))
data.pi    <-  cbind(y.pi, x.pi, union = TRUE)
end.period <- length(y.pi)
PI.EXP <- rep(NA,1,end.period+6)
for (t in seq(4,(end.period-39),1)){
    Y <- data.pi[seq(t,t+39,1),1]
    X <- cbind(as.matrix(data.pi[seq(t,t+39,1),c(2,3,4)]),
               rep(1,length(Y)))
    beta <- solve(t(X)%*%X)%*%t(X)%*%Y
    PI.EXP[t+48] <- c(as.matrix(data.pi[t+43,c(2,3,4)]),1)%*%beta
}
inflation.expectations <- tis(PI.EXP,
                              start = shiftQuarter(start(as.ts(x.pi)),-12),
                              tif = 'quarterly')

## Express FRBNY discount rate data on a 365-day basis
frbny.ffr.eff <- 100*((1+frbny.ffr/36000)^365 -1)

## NY Fed discount rate is used prior to 1965; thereafter, use the effective federal funds rate
interest <- mergeSeries(window(frbny.ffr.eff, end = c(1964,4)),
                        window(ffr, start = c(1965,1)))


##------------------------------------------------------------------------------##
## Output Data
##------------------------------------------------------------------------------##
data.out <- window(cbind(gdp.log,
                         inflation,
                         inflation.expectations,
                         oil.price.inflation,
                         import.price.inflation,
                         interest),
                   start = data.out.start)
write.table(data.out,file = 'LW_input_data.csv', sep = ',',
            col.names = TRUE, quote = FALSE, na = '.', row.names = FALSE)
