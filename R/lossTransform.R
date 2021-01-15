#' @title Linear Approximation of the object function
#' @name  lossTrans
#' @description Linear Approximation of the object function
#' @param obj the regression object from R output
#' @param n the sample size 
#' @import lars
#' @importFrom stats as.formula coef cor vcov
#' @rdname lossTransform
mod_lsa <- function(obj,n)
{
  
  if(class(obj)[1]=='coxph'|class(obj)=='fcrr') intercept <- 0

    
  Sigma <- vcov(obj)
  SI <- solve(Sigma)
  beta.ols <- coef(obj)
  l.fit <- lars.lsa(SI, beta.ols, intercept, n)

  t1 <- sort(l.fit$BIC, ind=T)
  t2 <- sort(l.fit$AIC, ind=T)
  beta <- l.fit$beta
  if(intercept) {
    beta0 <- l.fit$beta0+beta.ols[1]
    beta.bic <- c(beta0[t1$ix[1]],beta[t1$ix[1],])
    beta.aic <- c(beta0[t2$ix[1]],beta[t2$ix[1],])
  }
  else {
    beta0 <- l.fit$beta0
    beta.bic <- beta[t1$ix[1],]
    beta.aic <- beta[t2$ix[1],]
  }

  obj <- list(beta.ols=beta.ols, beta.bic=beta.bic,
              beta.aic = beta.aic)
  obj
}

grp_fun<-function(pmat, g) {  
  #### divide x into g groups
  px=length(pmat)
  ####specify how many groups
  k=floor(px/g)
  xgrp<-NULL
  for( m in 1:(g) ) 
    
  {cgrp<-data.frame(varname=pmat[(m*k-k+1): (m*k)] , grp=m)
  xgrp<-rbind(xgrp, cgrp)}
  remain<-px-k*g
  #cat('\nRemaining ',remain,' variables\n')
  if(px>k*g) {
    
    dgrp<-data.frame(varname=pmat[(k*g+1):px] , grp=sample(1:g, size=remain))
    xgrp<-rbind(dgrp, xgrp)
    
  }
  xgrp
}
crr.grp.boot<-function(data,index, myformula){
    set.seed(Sys.time())
    dat<-data[index, ]
    
   mycr<-fastcmprsk::fastCrr(myformula,data=dat, standardize = T)
    blsa<-mod_lsa(mycr, n=nrow(dat))
    myfirst<-blsa$beta.aic 
    names(myfirst)<-colnames(dat)[-c(1:2)]
    myfirst
  } 

###################################
## unified loss function transformation
lars.lsa <- function (Sigma0, b0, intercept,  n,
                      type = c("lasso", "lar"),
                      eps = .Machine$double.eps,max.steps) 
{
  type <- match.arg(type)
  TYPE <- switch(type, lasso = "LASSO", lar = "LAR")

  n1 <- dim(Sigma0)[1]
  
  ## handle intercept
  if (intercept) {
    a11 <- Sigma0[1,1]
    a12 <- Sigma0[2:n1,1]
    a22 <- Sigma0[2:n1,2:n1]
    Sigma <- a22-outer(a12,a12)/a11
    b <- b0[2:n1]
    beta0 <- crossprod(a12,b)/a11
  }
  else {
    Sigma <- Sigma0
    b <- b0
  }

  Sigma <- diag(abs(b))%*%Sigma%*%diag(abs(b))
  b <- sign(b)
  
  nm <- dim(Sigma)
  m <- nm[2]
  im <- inactive <- seq(m)
 
  Cvec <- drop(t(b)%*%Sigma)
  ssy <- sum(Cvec*b)
  if (missing(max.steps)) 
    max.steps <- 8 * m
  beta <- matrix(0, max.steps + 1, m)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ignores <- NULL

  while ((k < max.steps) & (length(active) < m)) {
    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if (!any(drops)) {
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
        R <- updateR(Sigma[inew, inew], R, drop(Sigma[inew, active]),
                     Gram = TRUE,eps=eps)
        if(attr(R, "rank") == length(active)) {
          ##singularity; back out
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action,  - inew)
        }
        else {
          if(first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
        }
      }
    }
    else action <- -dropid
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1
    
    if (length(active) >= m) {
      gamhat <- Cmax/A      
    }
    else {        
      a <- drop(w %*% Sigma[active, -c(active,ignores), drop = FALSE])
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }
    if (type == "lasso") {
      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      # cat('zmin ',zmin, ' gamhat ',gamhat,'\n') 
      if (zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
    }
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
    
    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w   
    Gamrat <- c(Gamrat, gamhat/(Cmax/A))

    arc.length <- c(arc.length, gamhat)
    if (type == "lasso" && any(drops)) {
      dropid <- seq(drops)[drops]
      for (id in rev(dropid)) {
        R <- downdateR(R,id)
      }
      dropid <- active[drops]
      beta[k + 1, dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
   
    actions[[k]] <- action
    inactive <- im[-c(active)]
  }
  beta <- beta[seq(k + 1), ]
 
  dff <- b-t(beta)

  RSS <- diag(t(dff)%*%Sigma%*%dff)

  if(intercept)
    beta <- t(abs(b0[2:n1])*t(beta))
  else
    beta <- t(abs(b0)*t(beta))
  
  if (intercept) {
    beta0 <- beta0-drop(t(a12)%*%t(beta))/as.vector(a11)
  }
  else {
    beta0 <- rep(0,k+1)
  }
  dof <- apply(abs(beta)>eps,1,sum)
  BIC <- RSS+log(n)*dof
  AIC <- RSS+2*dof
  
  object <- list(AIC = AIC, BIC = BIC, 
                 beta = beta, beta0 = beta0)
  object
}

lsa_glmnet <- function(obj, X, alpha=0.5){
  intercept <- attr(obj$terms,'intercept')
  if(class(obj)[1]=='coxph'|class(obj)[1]=='fcrr') intercept <- 0

Sigma <- vcov(obj)
SI <- solve(Sigma)
b0<-beta.ols <- coef(obj)

n1<-nrow(SI)
if (intercept) { 
  a11 <- SI[1,1]
  a12 <- SI[2:n1,1]
  a22 <- SI[2:n1,2:n1]
  SI.cond <- a22-outer(a12,a12)/a11
  b <- b0[2:n1]
  beta0 <- crossprod(a12,b)/a11
} else { SI.cond <-SI; b<-b0 }
X=as.matrix(X)
xtxInv=solve(t(X)%*%X)
xsand=X%*%xtxInv
Wmat=xsand%*%SI.cond%*%t(xsand)
r<-eigen(Wmat)
eignv<-r$vectors
lamgen<-r$values
lamgen[lamgen<=0]<-1e-6
w<-eignv%*%diag(lamgen^.5)%*%solve(eignv)

##adjust weight vi
vi=w #nrow(X)*
##new response variable
yita=vi%*%X%*%as.vector(b)
Xvi=vi%*%X
nd=nrow(X)
cvgnet<-glmnet::cv.glmnet(y=yita, x = Xvi, intercept=F,alpha=alpha, standardize=F, parallel = T)

nbeta=coef(cvgnet, s=cvgnet$lambda.1se) 
nbeta<-drop(nbeta)[-1]
nbeta

}





