##----------------------------------------------------------------------------------##
## File:        run.lw.R
##
## Description: This the main file for LW, which runs the three-stage LW estimation
##              for each economy and saves output.
##----------------------------------------------------------------------------------##
rm(list=ls())

##------------------------------------------------------------------------------##
## Load required packages and source all programs to be used in HLW estimation.
##------------------------------------------------------------------------------##
if (!require("tis")) {install.packages("tis"); library("tis")} ## Time series package
if (!require("nloptr")) {install.packages("nloptr"); library("nloptr")} ## Optimization

## Source all R programs
source("calculate.covariance.R")
source("format.output.R")
source("kalman.log.likelihood.R")
source("log.likelihood.wrapper.R")
source("utilities.R")
source("kalman.standard.errors.R")
source("kalman.states.R")
source("kalman.states.wrapper.R")
source("median.unbiased.estimator.stage1.R")
source("median.unbiased.estimator.stage2.R")
source("rstar.stage1.R")
source("rstar.stage2.R")
source("rstar.stage3.R")
source("unpack.parameters.stage1.R")
source("unpack.parameters.stage2.R")
source("unpack.parameters.stage3.R")


##------------------------------------------------------------------------------##
## Define variables
##------------------------------------------------------------------------------##
## Upper bound on a_3 parameter (slope of the IS curve)
a.r.constraint.stage2 <- -0.0025
a.r.constraint.stage3 <- -0.0025 

## Lower bound on b_2 parameter (slope of the Phillips curve)
b.y.constraint.stage1 <- NA
b.y.constraint.stage2 <- 0.025
b.y.constraint.stage3 <- 0.025

## Set the start and end dates of the estimation sample as well as the data start date (format is c(year,quarter))
data.start   <- c(1948,1)
sample.start <- c(1961,1) 
sample.end   <- c(2019,2)

## The estimation process uses data beginning 8 quarters prior to the sample start
est.data.start    <- shiftQuarter(sample.start,-8)

## Set start index for y for initialization of state vector
g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(est.data.start,'quarterly')

## Set column names for CSV output
output.col.names <- c("Date","rstar (filtered)","g (filtered)","z (filtered)","output gap (filtered)","","All results are output from the Stage 3 model.",rep("",12),"Standard Errors","Date","y*","r*","g","","rrgap","","Date","rstar (smoothed)","g (smoothed)","z (smoothed)","output gap (smoothed)")

## Set number of iterations for Monte Carlo standard error procedure
niter <- 5000

## Because the MC standard error procedure is time consuming, we include a run switch
## Set run.se to TRUE to run the procedure
run.se <- TRUE


##------------------------------------------------------------------------------##
## Read in data and compute inflation expectations series
##------------------------------------------------------------------------------##
data <- read.table("LW_input_data.csv",
                   sep=',',header=TRUE,stringsAsFactors=FALSE, na.strings=".")

## Get series beginning in data.start
log.output.all                      <- tis(data$gdp.log, start=data.start, tif='quarterly')
inflation.all                       <- tis(data$inflation, start=data.start, tif='quarterly')
relative.oil.price.inflation.all    <- tis(data$oil, start=data.start,
                                           tif='quarterly') - inflation.all
relative.import.price.inflation.all <- tis(data$import, start=data.start,
                                           tif='quarterly') - inflation.all
nominal.interest.rate.all           <- tis(data$interest, start=data.start, tif='quarterly')

inflation.expectations.all          <- tis(data$inflation.expectations, start=data.start, tif='quarterly')

## Get data in vector form beginning at est.data.start (set above)
log.output                      <- as.numeric(window(log.output.all, start=est.data.start))
inflation                       <- as.numeric(window(inflation.all, start=est.data.start))
relative.oil.price.inflation    <- as.numeric(window(relative.oil.price.inflation.all,
                                                     start=est.data.start))
relative.import.price.inflation <- as.numeric(window(relative.import.price.inflation.all,
                                                     start=est.data.start))
nominal.interest.rate           <- as.numeric(window(nominal.interest.rate.all, start=est.data.start))

inflation.expectations          <- as.numeric(window(inflation.expectations.all,
                                                     start=est.data.start))
real.interest.rate              <- nominal.interest.rate - inflation.expectations


##------------------------------------------------------------------------------##
## Run estimation
##------------------------------------------------------------------------------##

## Running the stage 1 model
out.stage1 <- rstar.stage1(log.output,
                           inflation,
                           relative.oil.price.inflation,
                           relative.import.price.inflation,
                           b.y.constraint=b.y.constraint.stage1)

## Median unbiased estimate of lambda_g
lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)

## Running the stage 2 model
out.stage2 <- rstar.stage2(log.output,
                           inflation,
                           relative.oil.price.inflation,
                           relative.import.price.inflation,
                           real.interest.rate,                         
                           lambda.g,
                           a.r.constraint=a.r.constraint.stage2,
                           b.y.constraint=b.y.constraint.stage2)

## Median unbiased estimate of lambda_z
lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)

## Running the stage 3 model
out.stage3 <- rstar.stage3(log.output,
                           inflation,
                           relative.oil.price.inflation,
                           relative.import.price.inflation,
                           real.interest.rate,
                           lambda.g,
                           lambda.z,
                           a.r.constraint=a.r.constraint.stage3,
                           b.y.constraint=b.y.constraint.stage3,
                           run.se)


##------------------------------------------------------------------------------##
## Save output
##------------------------------------------------------------------------------##

## One-sided (filtered) estimates
one.sided.est <- cbind(out.stage3$rstar.filtered,
                       out.stage3$trend.filtered,
                       out.stage3$z.filtered,
                       out.stage3$output.gap.filtered)

## Two-sided (smoothed) estimates
two.sided.est <- cbind(out.stage3$rstar.smoothed,
                       out.stage3$trend.smoothed,
                       out.stage3$z.smoothed,
                       out.stage3$output.gap.smoothed)

## Save output to CSV
output.us <- format.output(out.stage3, one.sided.est, two.sided.est, real.interest.rate, sample.start, sample.end, run.se = run.se)
write.table(output.us, 'LW_output.csv', col.names = output.col.names, quote=FALSE, row.names=FALSE, sep = ',', na = '')
