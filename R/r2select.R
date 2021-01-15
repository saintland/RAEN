#'@title Variable Selection with the candidate pool
#' @rdname r2select
#' @name r2select
#' @description Perform variable selection with pooled candidates
#' @param x.tr the predictor matrix
#' @param y.tr the time and status object for survival 
#' @param B times of bootstrap
#' @param weight variable weight 
#' @param prob variable selection probability
#' @param  parallel Logical TRUE or FALSE. Whether to use multithread computing, which can save consideratable amount of time for high dimensional data. Default is TRUE.
#' @param m the number of variables to be randomly included in the model in this step. Default is 8.  
#' @return the estimates of variables with B bootstraps, which is a dataframe with B rows and `ncol(x)` columns.  
#' @export 

r2select<-function(x.tr, y.tr, B, weight,prob, parallel=TRUE,m=8){ 
  if(parallel) `%mydo%` <- `%dopar%` else `%mydo%` <- `%do%`
  
  num.p<-ncol(x.tr)
  #m=ceiling(num.p*fr)
  nsize=nrow(x.tr)
  ## Bootstrap samples + random feature selection to assign weights
  
  ###parallel B boostrap
  parb<-grp_fun(1:B, g=20)
 j<-vector()
  par2<-foreach::foreach(j = c(1:20), .packages=c('fastcmprsk'),.verbose=TRUE) %mydo%{  
    pargr<-parb$varname[parb$grp==j]
    optbeta1 = rep(0,num.p)
    # indexmat1 = matrix(0,B,m)
    # somebetamat1 = matrix(0,B,m)
    denominator1 = rep(1e-8,num.p)
    beta.mat = matrix(0,B,num.p)
    for (b in pargr) {
      set.seed(Sys.time())
      samp.y = sample(1:nsize, nsize, replace=TRUE)
      
      samp.x = sample(1:num.p, m, replace=FALSE, prob = prob)
      #cat('\n',b,":",samp.x,'\n')
      
      ###check rlasso code
      newx.tr   = scale(x.tr[samp.y, samp.x], FALSE, weight[samp.x])
      #newx.tr   = x.tr[samp.y, samp.x]
      newy.tr   = y.tr[samp.y,]
      
      subdat<-data.frame(newy.tr,newx.tr)
      
      crm<-as.formula(paste('Crisk(time, cens)~', paste(names(subdat)[-c(1:2)],collapse = "+")) )
      mycr<-fastcmprsk::fastCrr(crm,data=subdat,standardize = T)
      
        mysecond<-lsa_glmnet(mycr, X=newx.tr)
            
      fbeta<-mysecond/weight[samp.x]
      optbeta1[samp.x] = optbeta1[samp.x] + fbeta
      denominator1[samp.x] = denominator1[samp.x] + rep(1,m)
      # indexmat1[b,] = sort(samp.x)
      # somebetamat1[b,] = fbeta[order(samp.x)]
      # beta.mat[b,sort(samp.x)]=fbeta[order(samp.x)]
      beta.mat[b,samp.x]=fbeta
      beta.mat[b,-samp.x]<-NA
      avebeta1 = optbeta1/denominator1
    }
    
    beta.mat[pargr,]
  }
  bmat<-do.call(rbind, par2)
  bmat
}

