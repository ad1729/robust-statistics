library(matrixcalc)
library(ccaPP)

rhoCauchy = function(x,b=0.306){
  # x is a univariate sample
  b2 = b*2
  rho = (x/b2)^2
  rho/(rho + 1)
} 

scale1StepM = function(x,rhofunction,precScale) {
  # Computes the first step of an algorithm for
  # a scale M-estimator using the given rho function. 
  # The scatter is computed relative to zero.
  # If you change rhofunction you must also change b.
  #
  x = x[!is.na(x)] # we always take out NAs
  n = length(x)
  if(n == 0) { return(0.0)
  } else {
    # sigma0 = 1.4826*median(abs(x))
    sigma0 = 1.4826*fastMedian(abs(x)) #=median deviation from zero * consistency factor voor normale verdeling
    if(sigma0 < precScale) { return(0.0)
    } else {
      rho = rhofunction(x/sigma0)
      return(sigma0 * sqrt(sum(rho)*2/n)) #deze twee is de "delta" van E_F(rho(x)) want E_F(rho(x))=1/2 als F std normaal
    }
  }
}

fastSplitSample = function(x){
  # Centers sample by median and deviations in 2 equal halves.
  # Assumes that NAs have already been removed.
  # This function has time complexity O(n).
  #
  med = fastMedian(x) # O(n) time
  x = x - med # centering
  n = length(x)  
  h = n %/% 2   # = integer part of n/2
  xa = x[x > 0] # length(xa) <= h
  xb = x[x < 0] # length(xa) <= h 
  xa = c(rep(0,(n - h - length(xa))),xa)
  xb = c(rep(0,(n - h - length(xb))),abs(xb)) #abs!!! 
  return(list(xa=xa,xb=xb,med=med))
}

compScales = function(x,robScale,rhofunction,rmZeroes,maxRatio,precScale){
  # assumes that x is an array of numbers
  x = x[!is.na(x)] # we always take out NAs
  # temp = splitSample(x)  
  temp = fastSplitSample(x)
  xa   = temp$xa
  xb   = temp$xb
  med  = temp$med
  sall = robScale((x - med),rhofunction=rhofunction,precScale=precScale)

  #if(sall < precScale) stop("exact fit situation")
  if(rmZeroes){ # reduces breakdown value but yields fewer implosions
    xa = xa[xa > 0]
    xb = xb[xb > 0]
  }
  sa = robScale(xa,rhofunction=rhofunction,precScale=precScale)
  sb = robScale(xb,rhofunction=rhofunction,precScale=precScale)
  if(!is.null(maxRatio)){
    if(maxRatio < 2) stop("maxRatio must be at least 2")
    sa = min(c(max(sa,sall/maxRatio,na.rm = TRUE),sall*maxRatio),na.rm = TRUE)
    sb = min(c(max(sb,sall/maxRatio,na.rm=TRUE),sall*maxRatio),na.rm = TRUE)
  } 
  return(list(sa=sa,sb=sb,med=med))
}

DO = function (data,robScale=scale1StepM,rhofunction=rhoCauchy,univ=FALSE,rmZeroes=FALSE,maxRatio=NULL,precScale=1e-12){
  #calculates the directional outlyingness of all points in a dataset
  # data in the form  (n x t x p) or (n x k x j x p) if multivariate
  # data in the form (n x t) or (n x k x j) if univariate
  # univ= TRUE if data is univariate
  
  if(univ){
    dims=dim(data)
    n=dims[1]
    t=dims[2]
    if (length(dims)==3){
      data=apply(data,c(1),as.vector)
      data=aperm(data,c(2,1))
      t=dims[2]*dims[3]
    }
    OUTPUT=apply(X =data ,MARGIN = 2,FUN= DO_univ,rhofunction=rhofunction,robScale=robScale,rmZeroes=rmZeroes,maxRatio=maxRatio,precScale=precScale)
    dim(OUTPUT)=c(n,2,t)
    don=OUTPUT[,1,]
    dod=OUTPUT[,2,]
    tempdo=don/dod
    f=function(x){#function to select scale where max is reached, if this scale exists
      out=0
      options(warn=-1)
      if(is.finite(max(x,na.rm=TRUE))){out=min(which(x==max(x,na.rm=TRUE)),na.rm = TRUE)}
      options(warn=0)
      return(out)
      }
    index=apply(tempdo,2, FUN = f)
    indices=cbind(index,1:t)
    lowerbound=dod[indices]
    lowerBound=exp((median(log(lowerbound[lowerbound>0]),na.rm = TRUE)-2*mad(log(lowerbound[lowerbound>0]),na.rm = TRUE)))
    
    dod=pmax(dod,lowerBound,na.rm=TRUE)
    DO_out=don/dod
    
    if(length(dims)==3){
      dim(DO_out)=dims[1:3] #keeps the data in the form n*k*j
      }
    return(DO_out)
  }
  else{
    dims=dim(data)
    n=dims[1]
    t=dims[2]
    p=dims[length(dims)]
    
    if (length(dims)==4){
      #when data in surface form -> unfold de pixels
      data=apply(data,c(1,4),as.vector)
      data=aperm(data,c(2,1,3)) #keeps dimensions in order n*t*p
      t=dims[2]*dims[3]
    }
    OUTPUT=apply(X = data,MARGIN = 2,FUN = ppsuit,robScale=robScale,rhofunction=rhofunction,lb=NULL,rmZeroes=rmZeroes,maxRatio=maxRatio,precScale=precScale)
    don=sapply(OUTPUT,"[[",1)
    dod=sapply(OUTPUT,"[[",2)
    
    lowerbound=as.vector(dod)
    lowerBound=exp((median(log(lowerbound[lowerbound>0]),na.rm = TRUE)-2*mad(log(lowerbound[lowerbound>0]),na.rm = TRUE)))
    
    zerodods=apply(dod,2,function(x) any(x==0,na.rm=TRUE)) #indices of those pixels containing zero dods => calculation has to be redone
    zerodods=which(zerodods==TRUE)
    
    dod=pmax(dod,lowerBound,na.rm=TRUE)
    DO_temp=don/dod  #t x n  this DO values are ok except for the pixels with 0 in the dod= zerodods pixels
    
    if (length(zerodods)!=0){
      data_zerodods=data[,zerodods,]
      if(length(zerodods)==1){dim(data_zerodods)=c(n,1,p)}
      OUTPUT_zerodods=apply(X = data_zerodods,MARGIN = 2,FUN = ppsuit,robScale=robScale,rhofunction=rhofunction,lb=lowerBound,rmZeroes=rmZeroes,maxRatio=maxRatio,precScale=precScale)
      
      don_zerodods=sapply(OUTPUT_zerodods,"[[",1)
      dod_zerodods=sapply(OUTPUT_zerodods,"[[",2)
      
      DO_zerodods=don_zerodods/dod_zerodods
      DO_temp[,zerodods]=DO_zerodods
    }
    DO_out=DO_temp
  
    if(length(dims)==4){
      dim(DO_out)=dims[1:3] #giet in de vorm n*k*j
      }
    return(DO_out)
  }
} 

ppsuit = function(X,lb=NULL,rhofunction,robScale,rmZeroes,maxRatio,precScale){
  # X of dim n x p
  n=dim(X)[1]
  p=dim(X)[2]
  ndir=250*p
  if (prod(apply(X,1,identical,X[1,]),na.rm=TRUE)==1){ #if all rows identical: return 0
    return(list(ppdon=rep(0,n), ppdod=rep(0,n)))
    
  }
  else{
    A=pickdir(X,ndir=ndir) #chooses exactly nrdir. directions
    
    Y <- X %*% t(A)
    out_temp=apply(X = Y,MARGIN = 2,FUN= DO_univ,rhofunction=rhofunction,robScale=robScale,rmZeroes=rmZeroes,maxRatio=maxRatio,precScale=precScale)
    
    dim(out_temp)=c(n,2,ndir)
    don_temp=out_temp[,1,] #n x ndir
    dod_temp=out_temp[,2,] #n x ndir
    
    if(!is.null(lb)){
      
      dod_temp=pmax(dod_temp,lb,na.rm = TRUE) # here we replace the dod by lb if dod<lb. We do this for all directions
    }
    tempdo=don_temp/dod_temp
    options(warn=-1)
    indexmax=unlist(apply(tempdo,1, function(x) min(which(x==max(x,na.rm=TRUE))))) #minimum takes first index in case of multiple
    options(warn=0)
    indix=cbind(1:n,indexmax) # indix[,j] gives coordinate (columnnr in dod) of the direction in which observation j has maximal DO 
    
    don_temp=don_temp[indix] # n x 1
    dod_temp=dod_temp[indix] # n x 1
    
    return(list(ppdon=don_temp, ppdod=dod_temp)) #don en dod: dim n x ndir
  }
}

DO_univ=function(x,robScale,rhofunction,rmZeroes,maxRatio,precScale){
  #x has to be univariate
  # calculates the directional outlyingness of all points in a univariate sample x
 
   y=array(NA,dim=c(length(x),2))
  
   temp=compScales(x,robScale=robScale,rhofunction=rhofunction,rmZeroes=rmZeroes,maxRatio=maxRatio,precScale=precScale)
    sa=temp$sa
    sb=temp$sb
    med=temp$med
   
  inda=which(x>med)
  indb=which(x<med)
  indm=which(x==med)

  y[inda,1]=x[inda]-med
  y[inda,2]=sa
  y[indb,1]=med-x[indb]
  y[indb,2]=sb
  y[indm,1]=0
  y[indm,2]=0 
  
  return(y)
  #return:column with numerators and denominators of DO
}

pickdir=function(X,ndir){
  # picks ndir normalized directions
  # dim X= n x p
  n=dim(X)[1]
  p=dim(X)[2]
  i=0 #counter for number of found directions
  j=0 # seed setter
  B=array(0.0,dim=c(ndir,p))
  while(i<ndir){
    j=j+1
    set.seed(j)
    indices=sample(n,p)
    Atemp=X[indices,]
    if(is.non.singular.matrix(Atemp,tol=1e-8)){
      nextdir=solve(Atemp)%*%matrix(rep(1,p)) # calculates the direction perpendicularto the hyperplane defined by the points in X[indices,]
      if (norm(as.matrix(nextdir),type="f")>1e-12){
        nextdir=nextdir/norm(as.matrix(nextdir),type="f")
        if(i==0){
          B[i+1,]=nextdir #rows of B are the random directions
          i=i+1
        }
        else{
          checkunique=B[1:i,]-matrix(rep(nextdir,i),byrow = TRUE,ncol = p) #check if we do not have this direction already
          mini=apply(checkunique,MARGIN=1,FUN = function(x) norm(as.matrix(x)))
          if(min(mini)>1e-12){
            B[i+1,]=nextdir #rows of B are the random directions
            i=i+1
          }
        }
      }
    }
  }
  return(B)
}

