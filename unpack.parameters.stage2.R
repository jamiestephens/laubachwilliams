##------------------------------------------------------------------------------##
## File:        unpack.parameters.stage2.R
##
## Description: This file generates coefficient matrices for the stage 2
##              state-space model for the given parameter vector.
##
## Stage 2 parameter vector: [a_1, a_2, a_3, a_4, a_5, b_1, b_2, b_3, b_4, b_5, sigma_y~, sigma_pi, sigma_y*]  
##------------------------------------------------------------------------------##
unpack.parameters.stage2 <- function(parameters, y.data, x.data, lambda.g, xi.00=NA, P.00=NA) {
  n.state.vars <- 4
    
  A         <- matrix(0, 2, 10)
  A[1, 1:2] <- parameters[1:2] ## a_1, a_2
  A[1, 3:4] <- parameters[3]/2 ## a_3
  A[1, 10 ] <- parameters[4]   ## a_4
  A[2, 1  ] <- parameters[8]   ## b_3
  A[2, 5:6] <- parameters[6:7] ## b_1, b_2
  A[2, 7  ] <- 1 - parameters[6] - parameters[7] ## b_1,b_2
  A[2, 8:9] <- parameters[9:10] ##b_4,b_5
  A         <- t(A)
  
  H         <- matrix(0, 2, 4)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2] ## a_1,a_2
  H[1, 4  ] <- parameters[5]    ## a_5
  H[2, 2]   <- -parameters[8]   ## b_3
  H         <- t(H)

  R         <- diag(c(parameters[11]^2, parameters[12]^2)) ## sigma_1,sigma_2
  Q         <- matrix(0, 4, 4)
  Q[1, 1]   <- parameters[13]^2              ## sigma_4
  Q[4, 4]   <- (lambda.g * parameters[13])^2 ## sigma_4

  F <- matrix(0, 4, 4)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- 1

  cons <- matrix(0, n.state.vars, 1)  
  
  ## Starting values for xi.00 and P.00
  if (any(is.na(xi.00))) {
      x  <- rbind(t(H), t(H) %*% F, t(H) %*% F %*% F, t(H) %*% F %*% F %*% F, t(H) %*% F %*% F %*% F %*% F)
      om <- matrix(0, 10, 10)
      om[5:6, 3:4] <- t(H) %*% F %*% Q %*% H
      om[7:8, 3:4] <- t(H) %*% F %*% F %*% Q %*% H
      om[9:10,3:4] <- t(H) %*% F %*% F %*% F %*% Q %*% H
      om[7:8, 5:6] <- t(H) %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om[9:10,5:6] <- t(H) %*% F %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om[9:10,7:8] <- t(H) %*% F %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H
      om           <- om + t(om)
      om[1:2, 1:2] <- R
      om[3:4, 3:4] <- t(H) %*% Q %*% H + R
      om[5:6, 5:6] <- t(H) %*% (F %*% Q %*% t(F) + Q) %*% H + R
      om[7:8, 7:8] <- t(H) %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H + R
      om[9:10, 9:10] <- t(H) %*% (F %*% F %*% F %*% Q %*% t(F) %*% t(F) %*% t(F) + F %*% F %*% Q %*% t(F) %*% t(F)
                                  + F %*% Q %*% t(F) + Q) %*% H + R
      p1 <- t(x) %*% solve(om, x)
      yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,], y.data[5,])
      tmp <- c(t(A) %*% x.data[1,],
               t(A) %*% x.data[2,],
               t(A) %*% x.data[3,],
               t(A) %*% x.data[4,],
               t(A) %*% x.data[5,])
      xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
      tmp <- yy - tmp - x %*% xi.00
      P.00 <- solve(p1, (diag(nrow=4) * sum(tmp^2) / (length(yy) - n.state.vars)))
  }
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data))
}
