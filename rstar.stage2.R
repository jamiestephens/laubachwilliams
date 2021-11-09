##------------------------------------------------------------------------------##
## File:        rstar.stage2.R
##
## Description: This file runs the model in the second stage of the LW estimation.
##------------------------------------------------------------------------------##
rstar.stage2 <- function(log.output,
                         inflation,
                         relative.oil.price.inflation,
                         relative.import.price.inflation,                         
                         real.interest.rate,
                         lambda.g,
                         a.r.constraint=NA,
                         b.y.constraint=NA,
                         xi.00=NA, P.00=NA) {

  stage <- 2
  
  ## Data must start 8 quarters before estimation period
  T <- length(log.output) - 8

  ## Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4), c(rep(0,56),1:(T+4-56)), c(rep(0,142),1:(T+4-142)))
  y.og <- log.output[5:(T+8)]
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

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
                  relative.import.price.inflation[9:(T+8)],
                  rep(1,T))
  
  ## Starting values for the parameter vector
  initial.parameters <- c(b.is, -b.is[3], b.ph[1:2], b.ph[4:6], s.is, s.ph, 0.5)

  ## Set an upper and lower bound on the parameter vectors:
  ## The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  ## Set a lower bound for the Phillips curve slope (b_2) of b.y.constraint, if not NA
  if (!is.na(b.y.constraint)) {
      print(paste("Setting a lower bound of b_y >",as.character(b.y.constraint),"in Stage 2"))
      if (initial.parameters[8] < b.y.constraint) {
          initial.parameters[8] <- b.y.constraint
      }
      theta.lb[8] <- b.y.constraint
  }

  ## Set an upper bound for the IS curve slope (a_3) of a.r.constraint, if not NA
  if (!is.na(a.r.constraint)) {
      print(paste("Setting an upper bound of a_r <",as.character(a.r.constraint),"in Stage 2"))
      if (initial.parameters[3] > a.r.constraint) {
          initial.parameters[3] <- a.r.constraint
      }
      theta.ub[3] <- a.r.constraint      
  }

  ## Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  theta <- nloptr.out$solution

  if (nloptr.out$status==-1 | nloptr.out$status==5) {
      print("Look at the termination conditions for nloptr in Stage 2")
      stop(nloptr.out$message)
  }

  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)$ll.cum

  ## Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)

  ## Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tt[,4] * 4
  potential.smoothed  <- c(states$smoothed$xi.tT[1, 3:2], states$smoothed$xi.tT[,1])
  output.gap.smoothed <- 100 * log.output[7:(T+8)] - potential.smoothed

  ## Inputs for median.unbiased.estimator.stage2.R
  y <- output.gap.smoothed[3:length(output.gap.smoothed)]
  x <- cbind(output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
             output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
             (x.data[,3]+x.data[,4])/2,
             states$smoothed$xi.tT[,4],
             rep(1,T))

  ## One-sided (filtered) estimates  
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

  ## Save variables to return
  return.list <- list()
  return.list$y              <- y
  return.list$x              <- x
  return.list$theta          <- theta
  return.list$log.likelihood <- log.likelihood
  return.list$states         <- states
  return.list$xi.00          <- xi.00
  return.list$P.00           <- P.00
  return.list$trend.filtered      <- trend.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$trend.smoothed      <- trend.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return(return.list)
}
