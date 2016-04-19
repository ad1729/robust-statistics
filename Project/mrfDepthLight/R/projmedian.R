projmedian <- function(x, projDepths=NULL, options=NULL){

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  #Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  #check projDepths
  if (!is.null(projDepths)) {
    projDepths <- data.matrix(projDepths)
    n2 <- nrow(projDepths)
    p2 <- ncol(projDepths)
    if (n1 != n2) { stop("A different number of depths from the number of observations was specified.") }
    if (p2 != 1) { stop("Provided depths have to be a columnvector.") }
    if (sum(projDepths > 1) != 0) { stop("The user supplied depths must take values in ]0,1].") }
    if (sum(projDepths <= 0) != 0) { stop("The user supplied depths must take values in ]0,1].") }
  }
  #check options
  if (is.null(options)) { options <- list() }
  if (!is.list(options)) { stop("options must be a list") }


  if (is.null(projDepths)) {
    projDepths <- projdepth(x = x, options = options)$depthX
  }

  return( findCenterProj(x,projDepths) )

}


findCenterProj <- function(x,projDepths){
  n <- nrow(x)
  p <- ncol(x)
  center <- vector("list",3)

  #MaxDepth
  maxDepth <- max(projDepths)
  indMaxDepths <- which(projDepths == maxDepth)
  center[[1]] <- colMeans(matrix(x[indMaxDepths,],ncol = p))
  names(center)[1] <- "max"

  #Gravity
  cutOff <- median(projDepths)
  indGravDepths <- which(projDepths >= cutOff)
  center[[2]] <- colMeans(matrix(x[indGravDepths,],ncol=p))
  names(center)[2] <- "gravity"


  #Huber
  weights <- matrix(rep(1,n),ncol=1)
  cutOff <- sqrt(qchisq(0.95,df=p))
  indHuber <- which( projDepths <=(1/(1+cutOff)) )
  outlyingness <- 1/projDepths[indHuber] -1
  weights[indHuber] <- (cutOff / outlyingness)^2
  center[[3]] <- t( (t(x) %*% weights)/n )
  names(center)[3] <- "Huber"

  return(center)
}

