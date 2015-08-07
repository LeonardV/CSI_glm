QLRT.Wald.Score.test.statistics <- function(fit0, fit1, s2, Amat) {

  beta0 <- coef(fit0)
  beta1 <- coef(fit1)
  X <- model.matrix(fit1)
  y <- model.response(fit1$model)
  n <- nrow(X)
  p <- ncol(X)
  I0 <- solve(vcov(fit0))
  I1 <- solve(vcov(fit1))
  
  QLR  <- 2/s2*( logLik(fit1)[1] - logLik(fit0)[1] )
  Wald <- 1/s2*(t(Amat%*%c(beta1-beta0)) %*% solve(Amat%*%solve(I1)%*%t(Amat)) %*% (Amat%*%c(beta1-beta0)))
  
  G0 <- colSums(score.vec(fit0))
  G1 <- colSums(score.vec(fit1))
  Score <- 1/s2*(G1-G0)%*%solve(I0)%*%(G1-G0)      
    
  cat(" ...QLR:",formatC(QLR, width=7, digits=4, format="f"), 
      " ...Wald:",formatC(Wald, width=7, digits=4, format="f"), 
      " ...Score:",formatC(Score, width=7, digits=4, format="f"), "\n")
  
  TS <- list(QLR=QLR, Wald=Wald, Score=Score)
  
  TS
}


#adjusted estfun() function from the sandwhich package
score.vec <- function (object, ...) {
  mm <- model.matrix(object)
  #the working residuals, that is the residuals in the final iteration of the IWLS fit. 
  #Since cases with zero weights are omitted, their working residuals are NA.
  wresid <- as.vector(residuals(object, "working"))*weights(object, "working")
  dispersion <- if (object$family[1] %in% c("poisson","binomial", 
                                       "Negative Binomial")) { 1 
                                                               } 
  else { 
    sum(wresid^2, na.rm=TRUE)/sum(weights(object, "working"), na.rm=TRUE) 
  }
  score.vec <- wresid * mm/dispersion
  
  score.vec
}
