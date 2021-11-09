##------------------------------------------------------------------------------##
## File:        unpack.parameters.stage1.R
##
## Description: This file generates coefficient matrices for the stage 1
##              state-space model for the given parameter vector.
##
## Stage 1 parameter vector: [a_1, a_2, b_1, b_2, b_3, b_4, b_5, g, sigma_1, sigma_2, sigma_4]  
##------------------------------------------------------------------------------##
unpack.parameters.stage1 <- function(parameters, y.data, x.data, xi.00=NA, P.00=NA) {
  n.state.vars <- 3
    
  A         <- matrix(0, 7, 2)
  A[1:2, 1] <- parameters[1:2]
  A[1, 2]   <- parameters[5]
  A[3:4, 2] <- parameters[3:4]
  A[5, 2]   <- 1-sum(A[3:4, 2])
  A[6:7, 2] <- parameters[6:7]
  
  H         <- matrix(0, 3, 2)
  H[1, 1]   <- 1
  H[2:3, 1] <- -parameters[1:2]
  H[2, 2]   <- -parameters[5]

  R         <- diag(c(parameters[9]^2, parameters[10]^2))
  Q         <- matrix(0, 3, 3)
  Q[1, 1]   <- parameters[11]^2

  F <- matrix(0, 3, 3)
  F[1, 1] <- F[2, 1] <- F[3, 2] <- 1

  cons <- matrix(0, 3, 1)
  cons[1, 1] <- parameters[8]
  
  if (any(is.na(xi.00))) {
      xi.00 <- rep(0, n.state.vars)
      
      ##  Starting values for xi.00 and P.00
      x  <- rbind(t(H), t(H) %*% F, t(H) %*% F %*% F, t(H) %*% F %*% F %*% F)
      om <- matrix(0, 8, 8)
      om[5:6, 3:4] <- t(H) %*% F %*% Q %*% H
      om[7:8, 3:4] <- t(H) %*% F %*% F %*% Q %*% H
      om[7:8, 5:6] <- t(H) %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om           <- om + t(om)
      om[1:2, 1:2] <- R
      om[3:4, 3:4] <- t(H) %*% Q %*% H + R
      om[5:6, 5:6] <- t(H) %*% (F %*% Q %*% t(F) + Q) %*% H + R
      om[7:8, 7:8] <- t(H) %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H + R
      
      p1 <- t(x) %*% solve(om, x) ## x' * inv(om) * x
      yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,])
      tmp <- c(t(A) %*% x.data[1,],
               (t(A) %*% x.data[2,] + t(H) %*% cons),
               (t(A) %*% x.data[3,] + t(H) %*% cons + t(H) %*% F %*% cons),
               (t(A) %*% x.data[4,] + t(H) %*% (diag(n.state.vars)+F+F%*%F) %*% cons))
      xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
      tmp <- yy - tmp - x %*% xi.00
      P.00 <- solve(p1, (diag(nrow=3) * sum(tmp^2) / 3))
  }
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data))
}
