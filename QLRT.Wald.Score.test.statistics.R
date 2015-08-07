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
  
  family <- fit1$family[1]
  link <- fit1$family[2]
  if(family=="gaussian") {
    ##Score test statistic
    #residuals under the eq. constrained model
    #I <- solve(1/(n-p) * as.numeric(t(res0)%*%res0) * solve(t(X)%*%X))
    s20 <- sum((y - X%*%beta0)^2)/(n-p)
    d0 <- t(1/s20 * t(X)%*%(y-X%*%beta0))
    d1 <- t(1/s2 * t(X)%*%(y-X%*%beta1))
    G0 <- d0
    G1 <- d1
    Score <- 1/s2*(G1-G0)%*%solve(I0)%*%t(G1-G0)      
  } 
  else if(family=="binomial") {
    if(link=="logit") {
      ##Score test statistic
      # The variance of a Bernouilli distribution is given by p(1-p)
      #pr = 1/(1+exp(-X%*%beta0))
      #V = array(0,dim=c(dim(X)[1],dim(X)[1]))
      #diag(V) = pr*(1-pr)
      #I = t(X)%*%V%*%X    
      G0 <- (y%*%X) - t(t(X)%*%(exp(X%*%beta0)/(1 + exp(X%*%beta0))))
      G1 <- (y%*%X) - t(t(X)%*%(exp(X%*%beta1)/(1 + exp(X%*%beta1))))
      Score <- 1/s2*(G1-G0)%*%solve(I0)%*%t(G1-G0)
    }
    else if(link=="probit") {
      ##Score test statistic
      #I <- matrix(0L, nrow=dim(X)[2], ncol=dim(X)[2])
      #for(i in 1:nrow(X)) {
      #  Z <- X[i,]
      #  Zb <- t(Z)%*%beta0
      #  I <- I + as.numeric((dnorm(Zb)^2 / ( pnorm(Zb)*(1-pnorm(Zb)) ))) * (Z%*%t(Z))
      #}
      phi0 <- pnorm(X%*%beta0); psi0 <- dnorm(X%*%beta0)
      phi1 <- pnorm(X%*%beta1); psi1 <- dnorm(X%*%beta1)
      G0 <- t(t(y*(psi0/phi0) - (1-y)*(psi0/(1-phi0)))%*%X)
      G1 <- t(t(y*(psi1/phi1) - (1-y)*(psi1/(1-phi1)))%*%X)
      Score <- 1/s2*t(G1-G0)%*%solve(I0)%*%(G1-G0)
    }
    else if(link=="log") {
      ##Score test statistic
      #I <- matrix(0L, nrow=dim(X)[2], ncol=dim(X)[2])
      #for(i in 1:nrow(X)) {
      #  Z <- X[i,]
      #  Zb <- t(Z)%*%beta0
      #  I <- I + as.numeric((dnorm(Zb)^2 / ( pnorm(Zb)*(1-pnorm(Zb)) ))) * (Z%*%t(Z))
      #}
      phi0 <- pnorm(X%*%beta0); psi0 <- dnorm(X%*%beta0)
      phi1 <- pnorm(X%*%beta1); psi1 <- dnorm(X%*%beta1)
      G0 <- t(t(y*(psi0/phi0) - (1-y)*(psi0/(1-phi0)))%*%X)
      G1 <- t(t(y*(psi1/phi1) - (1-y)*(psi1/(1-phi1)))%*%X)
      #G1 <- t(colSums(sum(y*(psi1/phi1) - (1-y)*(psi1/(1-phi1)))*X))
      Score <- 1/s2*t(G1-G0)%*%solve(I0)%*%(G1-G0)
    }
  }  
  else if(family=="poisson") {
    #Score test statistic
    #I <- matrix(0L, nrow=dim(X)[2], ncol=dim(X)[2])
    #for(i in 1:nrow(X)){
    #  bZ <- t(beta0)%*%(X[i,])
    #  Z <- X[i,]
    #  ZZ <- Z%*%t(Z)
    #  I <- I + ZZ*as.numeric(exp(bZ))
    #}
    G0 <- t(t(X)%*%(y-exp(X%*%beta0)))
    G1 <- t(t(X)%*%(y-exp(X%*%beta1)))
    Score <- 1/s2*(G1-G0)%*%solve(I0)%*%t(G1-G0)
  }
  
  cat(" ...QLR:",formatC(QLR, width=7, digits=4, format="f"), 
      " ...Wald:",formatC(Wald, width=7, digits=4, format="f"), 
      " ...Score:",formatC(Score, width=7, digits=4, format="f"), "\n")
  
  TS <- list(QLR=QLR, Wald=Wald, Score=Score)
  
  TS
  }
