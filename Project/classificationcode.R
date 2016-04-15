DOclassif=function(x,z,robScale=scale1StepM,rhofunction=rhoCauchy,rmZeroes=FALSE,maxRatio=NULL,precScale=1e-12){
# calculates the DO of z with respect to the points in x
# x of dimesions n1 x p, z of dimensions n2 x p
# also works for univariate datasets, i.e. x of dim n1 and z of dimension n2
if(is.null(dim(x))){
  n1=length(x)
  n2=length(z)
  temp=compScales(x=x,robScale = robScale,rhofunction = rhofunction,rmZeroes = rmZeroes,maxRatio = maxRatio,precScale = precScale)
  sa=temp$sa
  sb=temp$sb
  med=temp$med
  
  DOz=z
  DOz[DOz>=med]=(DOz[DOz>=med]-med)/sa
  DOz[DOz<med]=(med-DOz[DOz<med])/sb
  return(DOz)
  }
else{
n1=dim(x)[1]
p=dim(x)[2]
n2=dim(z)[1]
ndir=250*p
A=pickdir(x,ndir=ndir) #chooses exactly nrdir. directions, each direction is a row in A

Y <- x %*% t(A) #projected data, each column is the dataset projected on a direction

tempresult=apply(X = Y,MARGIN = 2,FUN = compScales,robScale = robScale,rhofunction = rhofunction,rmZeroes = rmZeroes,maxRatio = maxRatio,precScale = precScale)
medians=sapply(tempresult,"[[",3)
sas=sapply(tempresult,"[[",1)
sbs=sapply(tempresult,"[[",2)

Y2<-z%*%t(A) #columns contain projected dataset z
output=Y2
medians=t(array(medians,dim = c(ndir,n2)))
sas=t(array(sas,dim = c(ndir,n2)))
sbs=t(array(sbs,dim = c(ndir,n2)))
Y2[Y2>=medians]=(Y2[Y2>=medians]-medians[Y2>=medians])/sas[Y2>=medians]
Y2[Y2<medians]=(medians[Y2<medians]-Y2[Y2<medians])/sbs[Y2<medians]
DOz=apply(Y2, MARGIN = 1, function(x) max(x, na.rm=TRUE))
return(DOz)
}
}