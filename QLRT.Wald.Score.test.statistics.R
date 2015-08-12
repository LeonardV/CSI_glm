QLRT.Wald.Score.test.statistics <- function(fit0, fit1, s2, verbose=TRUE) {

  Amat <- fit0$Amat
  beta0 <- coef(fit0)
  beta1 <- coef(fit1)
  X <- model.matrix(fit1)
  y <- model.response(fit1$model)
  n <- nrow(X)
  p <- ncol(X)
  I0 <- solve(vcov(fit0))
  I1 <- solve(vcov(fit1))
  
  dispersion <- if (fit0$family[1] %in% c("poisson","binomial", "Negative Binomial")) { 
    1L 
  } 
  else { 
    s2
  }
  
  #Quasi likelihood ratio test statistics
  QLR  <- 2/dispersion*( logLik(fit1)[1] - logLik(fit0)[1] )
  
  #Wald test statistic
  Wald <- 1/dispersion*(t(Amat%*%c(beta1-beta0)) %*% solve(Amat%*%solve(I1)%*%t(Amat)) %*% (Amat%*%c(beta1-beta0)))
  #Wald <- 1/s2*(t(Amat%*%beta1) %*% solve(Amat%*%solve(I1)%*%t(Amat)) %*% (Amat%*%beta1))
  
  #Score test statistic
  G0 <- colSums(score.vec(fit0))
  G1 <- colSums(score.vec(fit1))
  Score <- 1/dispersion*(G1-G0)%*%solve(I0)%*%(G1-G0)      
  
  #Generalized distance test statistic
#  d1 <- t(beta2-beta0) %*% I2 %*% (beta2-beta0) - t(beta2-beta1) %*% I2 %*% (beta2-beta1)
#  d2 <- t(Amat%*%beta2-Amat%*%beta0) %*% (solve(Amat%*%solve(I2)%*%t(Amat))) %*% (Amat%*%beta2-Amat%*%beta0) - 
#    t(Amat%*%beta2-Amat%*%beta1) %*% (solve(Amat%*%solve(I2)%*%t(Amat))) %*% (Amat%*%beta2-Amat%*%beta1)
  
#  GD1 <- d1
#  GD2 <- d2
  
  #Pearson chi-square test statitic for grouped data
#  X2 <- t((y-X%*%beta0)^2 - (y-X%*%beta1)^2) %*% (1/(X%*%beta0)) #CHECKME!
  
  
  if(verbose) {
    cat(" ...QLR:",formatC(QLR, width=7, digits=4, format="f"), 
        " ...Wald:",formatC(Wald, width=7, digits=4, format="f"), 
        " ...Score:",formatC(Score, width=7, digits=4, format="f"), "\n")
#        " ...Generalized D1:",formatC(GD1, width=7, digits=4, format="f"),
#        " ...Generalized D2:",formatC(GD2, width=7, digits=4, format="f"), "\n")
  }
  
  TS <- list(QLR=QLR, Wald=Wald, Score=Score)#, GD1=GD1, GD2=GD2)
  
  TS
}


#adjusted estfun() function from the sandwhich package
score.vec <- function (object, ...) {
  mm <- model.matrix(object)
  #the working residuals, that is the residuals in the final iteration of the IWLS fit. 
  #Since cases with zero weights are omitted, their working residuals are NA.
  wresid <- as.vector(residuals(object, "working"))*weights(object, "working")
  dispersion <- if (object$family[1] %in% c("poisson","binomial", "Negative Binomial")) { 
    1 
  } 
  else { 
    sum(wresid^2, na.rm=TRUE)/ (sum(weights(object, "working"), na.rm=TRUE) - 0) #likelihood method
  }
  score.vec <- (wresid * mm) / dispersion
  
  score.vec
}
