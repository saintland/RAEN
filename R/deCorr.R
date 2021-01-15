#' @title De-correlating variables
#' @name deCorr
#' @description Divide the highly correlated variables into exclusive groups
#' @param x the predictor matrix
#' @param ngrp the number of blocks to separate variables 
#' @param rho the preset correlation threshold. Variables with correlation higher than rho will be separate into exclusive groups. Default is set to 0.7
#' @return a dataframe of variable names `varname` and the variable subgroup membership `grp` 
#' @export 
#' @rdname deCorr
deCorr<-function(x, rho=.7,ngrp=floor(15*ncol(x)/nrow(x))){

 nsize=nrow(x)
 psize=ncol(x)
  
   xcor<-abs(cor(x))

 sqcor<-apply(xcor, 1, function(x) sum(x[x<1 & x>0.5], na.rm=T) )

 xmax<-which.max(sqcor)
 xvec<-xcor[xmax,]
 rhot<-rho
 xclose<-xvec[xvec>=rhot]
 xrem<-which(xvec<rhot)
 xlen<-length(xclose)
 myclust<-list()
 myclust[[1]]<-xclose
 jj=1
 cat( '\nNo:',jj, 'cluster contains :', xlen,', remaining', length(xrem),'\n')
 while (xlen>0) {
   jj=jj+1
   rcor<-sqcor[names(sqcor) %in% names(xrem)]
   ymax<-which.max(rcor)
   yvec<-xcor[rownames(xcor) %in% names(xrem),][ymax,]
   yvec<-yvec[setdiff(names(yvec),names(xclose))]
   myclust[[jj]]<-yvec[yvec>=rhot]
   xclose<-c(xclose, yvec[yvec>=rhot])
   xlen<-length(yvec[yvec>=rhot])
   xrem<-yvec[yvec<rhot]
   cat('\nNo:',jj, 'cluster contains :', xlen,', remaining', length(xrem),'\n')
   if(xlen<5) break('Less than 5 rho>= threshold. De-Correlation stopped')
 }
clustlen<-sapply(myclust, length)
clustlen<-sort(clustlen, decreasing=T)
ngrp<-ifelse(ngrp<max(clustlen), clustlen[2], ngrp)

  bgrp<-sample(names(xrem), size=length(xrem), replace = F)
  base_grp<-grp_fun(bgrp, g=ngrp)
  base_grp$varname<-as.character(base_grp$varname)
  cor_break<-function(myclust){
  grpout<-NULL
  for (i in 1:jj) {
    grpk<-length(myclust[[i]])
    km<-sample(1:ngrp, size=grpk)
    
    grpm<-data.frame(varname=names(myclust[[i]]), grp=km)
    grpout<-rbind(grpout, grpm )
  }
  grpout}
  
  corg<-cor_break(myclust)
  ctable<-table(corg$grp)
  corg$varname<-as.character(corg$varname)
  fl<-1
  
  while(fl>0){
    for (i in names(ctable[ctable>1]))
    { vg<-corg$varname[corg$grp==i]
    
    srow<-xcor[rownames(xcor) %in% vg[1],, drop=F]
    xycor<-srow[,colnames(srow) %in% vg[-1]]
    if(any(xycor>rhot)) {
#cat('\n',vg,i,'\n');
break; corg<-cor_break(myclust); ctable<-table(corg$grp)}
    } 
    fl=0}
    
  ###combine the backgroup variable groups with correlated variable groups
  fgrp<-rbind(corg, base_grp)
  fgrp}