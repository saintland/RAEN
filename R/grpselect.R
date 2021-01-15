#' @title grpselect
#' @name grpselect
#' @description This is the split step, where variable in subgroups are selected
#' @param fgrp the variable group object from `deCorr``
#' @param x the predictor matrix
#' @param y a dataframe of time to event and event status. The primary outcome status is coded 1, the secondary outcome as 2, etc. The censored is coded as 0. 
#' @param B the number of bootstraps 
#' @param parallel whether to use multiple cores for parallel computing. Default is TRUE. 
#' @return a list of 
#' \itemize{
#' \item{fselect:}{ Names of the selected variables.}
#' \item{prob:}{ the generalized ridge variable importance.}
#' \item{weight:}{ the inverse of the ridge variable importance.}
#'} 
#' @export 
#' @rdname grpselect


grpselect<-function(fgrp, x,y,B=50, parallel=TRUE) {
  if(parallel) `%mydo%` <- `%dopar%` else  `%mydo%` <- `%do%`
  
  grp<-unique(fgrp$grp)
j<-vector()
  bt.res<-foreach::foreach(j = grp, .packages=c('fastcmprsk'), .verbose=TRUE) %mydo%{ 
     
    grpvar=fgrp$varname[fgrp$grp==j]
    nvar=length(grpvar)
    
    xformula<-paste(grpvar, collapse = '+')
    
    cat('\ngroup :',j,'nvar:', nvar, xformula, '\n')
    data=data.frame(y,x)
    subdat<-data[, c('time','cens', grpvar)]
    crm<-as.formula(paste('Crisk(time, cens)~', paste(names(subdat)[-c(1:2)],collapse = "+")) )
    bt.est<-boot::boot(data=subdat,statistic=crr.grp.boot, myformula= crm, R=B)
    coefmat<-bt.est$t
    dim(coefmat)
    colnames(coefmat)<-grpvar
    coefmat}
  
  fbt.grp<-do.call('cbind',bt.res)
  fbt.grp<-fbt.grp[,colnames(x)]
  
  
  
  small.cut<-1/nrow(x)
  fave<- apply(fbt.grp,2, function(x) sum(abs(x)>small.cut)/B )
  #### step1, set alpha, I_jk
  alpha.first=.95

  cat('\n', length(first.beta<-apply(fbt.grp[,fave>=alpha.first],2,mean)), 'predictors go to round two','\n')
  fselect<-colnames(x)[fave>=alpha.first]
  step1beta<-apply(fbt.grp[,fselect],2,mean)
  ####use beta from Round 1 as weight for Round 2
  x1mat<-as.matrix(x[,fselect])
  xtx<-t(x1mat)%*%x1mat
  bxtx<-t(step1beta)%*%xtx%*%step1beta+1e-8
  fprob=as.vector(xtx%*%step1beta)*step1beta/bxtx
  fprob<-abs(fprob)
  fweight<-1/(fprob+1e-8)
  
  list(weight=fweight, prob=fprob, fselect=fselect)
}