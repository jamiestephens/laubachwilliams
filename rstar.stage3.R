##------------------------------------------------------------------------------##
## File:        rstar.stage3.R
##
## Description: This file runs the model in the third stage of the LW estimation.
##------------------------------------------------------------------------------##
rstar.stage3 <- function(log.output,
                         inflation,
                         relative.oil.price.inflation,
                         relative.import.price.inflation,
                         real.interest.rate,
                         lambda.g,
                         lambda.z,
                         a.r.constraint=NA,
                         b.y.constraint=NA,
                         run.se=FALSE,
                         xi.00=NA, P.00=NA) {

  stage <- 3
  
  ## Data must start 8 quarters before estimation period
  T <- length(log.output) - 8

  ## Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4), c(rep(0,56),1:(T+4-56)), c(rep(0,142),1:(T+4-142)))
  y.og <- log.output[5:(T+8)]
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

  
  ## Initialization of state vector for Kalman filter using HP trend of log output
  ## Pulled into unpack.parameters.stage3.R
  b.pot <- solve(t(x.og) %*% x.og, t(x.og) %*% y.og)
  g.pot <- x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)
  xi.00.gpot <- c(100*g.pot[5:3],100*b.pot[2],100*b.pot[2],0,0)
                                               
  ## IS curve
  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)],
                (real.interest.rate[8:(T+7)] + real.interest.rate[7:(T+6)])/2,
                rep(1,T))
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is)
  r.is <- y.is - x.is %*% b.is
  s.is <- sqrt(sum(r.is^2) / (T-dim(x.is)[2]))

  ## Phillips curve
  y.ph <- inflation[9:(T+8)]
  x.ph <- cbind(inflation[8:(T+7)],
                (inflation[7:(T+6)]+inflation[6:(T+5)]+inflation[5:(T+4)])/3,
                (inflation[4:(T+3)]+inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/4,
                output.gap[4:(T+3)],
                relative.oil.price.inflation[8:(T+7)],
                relative.import.price.inflation[9:(T+8)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (T-dim(x.ph)[2]))
  
  y.data <- cbind(100 * log.output[9:(T+8)],
                  inflation[9:(T+8)])
  x.data <- cbind(100 * log.output[8:(T+7)],
                  100 * log.output[7:(T+6)],
                  real.interest.rate[8:(T+7)],
                  real.interest.rate[7:(T+6)],                  
                  inflation[8:(T+7)],
                  (inflation[7:(T+6)]+inflation[6:(T+5)]+inflation[5:(T+4)])/3,
                  (inflation[4:(T+3)]+inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/4,
                  relative.oil.price.inflation[8:(T+7)],
                  relative.import.price.inflation[9:(T+8)])
  
  initial.parameters <- c(b.is[1:3], b.ph[1:2], b.ph[4:6], 1, s.is, s.ph, 0.7)
  
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  ## Set a lower bound for the Phillips curve slope (b_2) of b.y.constraint, if not NA
  if (!is.na(b.y.constraint)) {
      print(paste("Setting a lower bound of b_y >",as.character(b.y.constraint),"in Stage 3"))
      if (initial.parameters[6] < b.y.constraint) {
          initial.parameters[6] <- b.y.constraint
      }
      theta.lb[6] <- b.y.constraint
  }

  ## Set an upper bound for the IS curve slope (a_3) of a.r.constraint, if not NA
  if (!is.na(a.r.constraint)) {
      print(paste("Setting an upper bound of a_r <",as.character(a.r.constraint),"in Stage 3"))
      if (initial.parameters[3] > a.r.constraint) {
          initial.parameters[3] <- a.r.constraint
      }
      theta.ub[3] <- a.r.constraint      
  }

  P.00.gpot <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00=xi.00.gpot)
    
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage,
                                                       lambda.g, lambda.z,
                                                       xi.00=NA, P.00=NA,
                                                       xi.00.gpot, P.00.gpot)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  theta <- nloptr.out$solution
  
  if (nloptr.out$status==-1 | nloptr.out$status==5) {
      print("Look at the termination conditions for nloptr in Stage 3")
      stop(nloptr.out$message)
  }

  ## Get xi.00 and P.00
  init.vals <- unpack.parameters.stage3(theta, y.data, x.data, lambda.g, lambda.z,
                                        xi.00=NA, P.00=NA, xi.00.gpot, P.00.gpot)
  xi.00 <- init.vals$xi.00
  P.00  <- init.vals$P.00
  print(paste("GLS switch value:",init.vals$gls))
  
  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

  ## If run.se = TRUE, compute standard errors for estimates of the states and report run time
  if (run.se) {
      ptm <- proc.time()
      se <- kalman.standard.errors(T, states, theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00, niter, a.r.constraint, b.y.constraint)
      print("Standard error procedure run time")
      print(proc.time() - ptm)
  }
  
  ## One-sided (filtered) estimates
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  z.filtered          <- states$filtered$xi.tt[,6]
  rstar.filtered      <- trend.filtered * theta[9] + z.filtered
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

  ## Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tT[,4] * 4
  z.smoothed          <- states$smoothed$xi.tT[,6]
  rstar.smoothed      <- trend.smoothed * theta[9] + z.smoothed
  potential.smoothed  <- states$smoothed$xi.tT[,1]/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)
   
  ## Save variables to return
  return.list <- list()
  return.list$rstar.filtered      <- rstar.filtered
  return.list$trend.filtered      <- trend.filtered
  return.list$z.filtered          <- z.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$rstar.smoothed      <- rstar.smoothed
  return.list$trend.smoothed      <- trend.smoothed
  return.list$z.smoothed          <- z.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return.list$theta               <- theta
  return.list$log.likelihood      <- log.likelihood
  return.list$states              <- states
  return.list$xi.00               <- xi.00
  return.list$P.00                <- P.00
  return.list$y.data              <- y.data
  return.list$initial.parameters  <- initial.parameters
  return.list$init.vals           <- init.vals
  if (run.se) { return.list$se    <- se }
  return(return.list)
}
