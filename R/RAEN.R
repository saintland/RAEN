#' @title Random Ensemble Variable Selection for Ultra-high Dimensional Data
#' @name  RAEN
#' @description Perform variable selection for ultra-high dimensional data
#' @param x the predictor matrix
#' @param y the time and status object for survival 
#' @param B times of bootstrap
#' @param ngrp the number of blocks to separate variables into. Default is 15*p/N, where p is the number of predictors and N is the sample size.
#' @param  parallel Logical TRUE or FALSE. Whether to use multithread computing, which can save consideratable amount of time for high dimensional data
#' @param ncore Number of cores used for parallel computing, if parallel=TRUE
#' @return a dataframe with the variable names and the regression coefficients
#' @rdname RAEN
#' @export
RAEN<-function(x,y,B,ngrp=floor(15*ncol(x)/nrow(x)), parallel=TRUE, ncore=2){
  
  if (parallel){
  myClust <- parallel::makeCluster(ncore)
  doParallel::registerDoParallel(myClust)
  }
  # Are we using multiple cores (parallel) or not
  if(parallel) `%mydo%` <- `%dopar%` else  `%mydo%` <- `%do%`
  nsize=nrow(x)
  psize=ncol(x)
  data=data.frame(y,x)
  colnames(data)[1:2]<-c('time', 'cens')
  ####find the assignment of variables
  fgrp<-deCorr(x, ngrp=ngrp)
  ####variable selection in groups
  f1p<-grpselect(fgrp, x,y,B=50)

fselect<-f1p$fselect
n.ratio<-length(fselect)/ncol(x)
num.first<-length(fselect)
fr= 8/num.first 
B=500

cat('\n','fraction is ', fr,'\n')
fweight<-f1p$weight
fprob<-f1p$prob
est.mat<-r2select(x.tr = data[, fselect], y.tr = data[,c('time', 'cens') ], B=B, fr=fr, weight = fweight, prob=fprob)

est.mat<-data.frame(est.mat)
colnames(est.mat)<-fselect
stopCluster(myClust)

small.cut<-1/nrow(x)
  beta.ave<-NULL
  beta.ave<-apply(est.mat, 2, function(x) { x<-x[!is.na(x)]
  denom<-sum(x!=0)
  if (denom==0|is.na(denom)) est=0 else est<-mean(x,na.rm=T)
  est})
names(beta.ave)<-fselect
beta.ave[abs(beta.ave)<=small.cut]<-0
beta.ave[is.na(beta.ave)]<-0
bout<- beta.ave[abs(beta.ave)>small.cut]


raen<-data.frame( id=names(bout), coef=bout)
#class(raen)<-'RAEN'
raen
}

#' @rdname RAEN
#' @method predict RAEN
#' @param obj the RAEN object containing the variable selection results
#' @param newdata the predictor matrix for prediction
#' @return the linear predictor of the outcome risk
#' @export

predict.RAEN<-function(obj,newdata,...) {
  inx<-newdata[,match(obj$id, colnames(newdata))]
  raen.lp<-as.matrix(inx)%*%as.numeric(obj$coef)
  raen.lp}

  