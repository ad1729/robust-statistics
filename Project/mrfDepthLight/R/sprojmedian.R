sprojmedian<-function(x, sprojDepths=NULL, options=NULL) {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  #Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n1 <- nrow(x)
  if (n1 > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  #check spDepths
  if (!is.null(sprojDepths)) {
    sprojDepths <- data.matrix(sprojDepths)
    n2 <- nrow(sprojDepths)
    p2 <- ncol(sprojDepths)
    if (n1 != n2) { stop("A different number of depths from the number of observations was specified.") }
    if (p2 != 1) { stop("Provided depths have to be a columnvector.") }
    if (sum(sprojDepths > 1) != 0) { stop("The user supplied depths must take values in ]0,1].") }
    if (sum(sprojDepths <= 0) != 0) { stop("The user supplied depths must take values in ]0,1].") }
  }
  #check options
  if (is.null(options)) { options <- list()  }
  if (!is.list(options)) { stop("options must be a list") }


  if (is.null(sprojDepths)) {
    sprojDepths <- sprojdepth(x = x, options = options)$depthX
  }

  return( findCenter(x,sprojDepths) )

}


findCenter <- function(x,sprojDepths) {
  p <- ncol(x)
  center <- vector("list",2)

  #MaxDepth
  maxDepth <- max(sprojDepths)
  indMaxDepths <- which(sprojDepths == maxDepth)
  center[[1]] <- colMeans(matrix(x[indMaxDepths,],ncol = p))
  names(center)[1] <- "max"

  #Gravity
  cutOff <- median(sprojDepths)
  indGravDepths <- which(sprojDepths >= cutOff)
  center[[2]] <- colMeans(matrix(x[indGravDepths,],ncol = p))
  names(center)[2] <- "gravity"

  return(center)
}

