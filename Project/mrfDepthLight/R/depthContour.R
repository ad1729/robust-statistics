depthContour <- function(x, alpha=NULL, type="hdepth", directions=NULL, options=NULL){

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  #Check the x data.
  x <- data.matrix(x)
  if (sum(is.nan(x) != 0)) { stop("Missing values in matrix x are not allowed.") }
  if (!is.numeric(x)) { stop("The input argument x must be a numeric data matrix.") }
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2) { stop('x should be at least two dimensional.') }
  if (n < (p + 1))			stop("The should a least be 3 observations.")

  #check alpha
	if (!is.numeric(alpha)) { stop("alpha should be a numeric vector.") }
	if (sum(sort(alpha) == alpha) != length(alpha)) { stop('The vector alpha should be sorted in ascending order.') }

  #Check type
  Indtype <- match(type,c("hdepth","projdepth","sprojdepth"))[1]
  if (is.na(Indtype)) { stop("type should be one of hdepth, projdepth or sprojdepth.") }
  if (min(alpha) < 1/n && Indtype == 1) {  stop("alpha should be at least 1/n for halfspace depth.") }
  if (max(alpha) > 0.5 && Indtype == 1) {	stop("alpha should be < 0.5.") }
  if (min(alpha) <= 0 || max(alpha) >= 1  ) { stop("alpha should lie in the open interval (0,1).") }

  #Check directions
  if (is.null(directions)) {
      set.seed(123)
      NDirections <- 250*p
      directions <- matrix(rnorm(n = NDirections * p),ncol = p)
      directions <- directions / sqrt( rowSums(directions^2) )
  } else{
    directions <- data.matrix(directions)
    if (sum(is.nan(directions) != 0)) { stop("Missing values in matrix directions are not allowed.") }
    if (!is.numeric(directions)) { stop("The input argument directions must be a numeric data matrix.") }
    pDir <- ncol(directions)
    if (pDir != p) { stop("The dimension of the directions vector should be the same as x.")}
  }
  directions <- directions / rowSums(directions^2)

  #check options
  if (is.null(options)) { options <- list(approx = TRUE) }
  if (!is.list(options)) { stop("options must be a list") }
  if ( "maxIter"  %in% names(options) ) {
    maxIter <- options[["maxIter"]]
    if (is.numeric(maxIter)) {
      if (maxIter < 1) { stop("Option maxIter must be a strictly positive integer.") }
    }
    if (!is.numeric(maxIter)) {
      stop("Option maxIter must be a strictly positive integer.")
    }
  } else {
    options$maxIter <- 100
  }


  #Check data for exact fits
  tol <- 1e-7
  ScaledX <- scale(x)
  temp <- attributes(ScaledX)
  ColumnSd <- temp[["scaled:scale"]]
  if (sum(ColumnSd <= 1e-14) > 0) {
    warning("One of the variables has deviation zero. Check the data matrix x.")
    return(list(depth = NULL,
                singularSubsets = NULL,
                dimension = sum(ColumnSd > 1e-14),
                hyperplane = as.numeric(ColumnSd <= 1e-14),
                type = "hdepth"
    )
    )
  }
  w1 <- try(svd(ScaledX/sqrt(n - 1)),silent = TRUE)
  if (!is.list(w1)) {
	  warning("The singular-value decomposition of the data matrix x could not be computed.")
	  Result <- list(vertices = NULL,
                 depth = NULL,
                 count = NULL,
                 empty = NULL,
                 dithered = NULL,
                 dimensions = NULL,
                 hyperplane = NULL,
	               type = "hdepthContour")
		return(Result)
	}
	if (min(w1$d) < tol) {
		warning("An exact fit was found. Check output for more details.")
		Result <- list(vertices = NULL,
		             depth = NULL,
		             count = NULL,
		             empty = NULL,
		             dithered = NULL,
		             dimensions = sum(w1$d > tol),
                 hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
		             type = "hdepthContour")
    return(Result)
	}

  if (type == "hdepth") {
    if (p == 2) { Result <- CalcBivHContour(x, alpha) }
    if (p > 2) { Result <- CalcHContour(x, alpha, directions, options) }
  }
  else if (type == "projdepth" || type == "sprojdepth") {
    Result <- CalcProjTypeContour(x, alpha, type, directions, options)
  }
  else{
    stop("type not recognized.")
  }

  return(Result)
}

CalcBivHContour <- function(x, alpha) {

  n <- nrow(x)

  Result <- vector("list",length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
      TResult <- .Fortran("iso2hdw",
                          as.double(x[,1,drop = TRUE]),  #1 First variable of the data set.
                          as.double(x[,2,drop = TRUE]),  #2 Second variable of the data set.
                          as.integer(n),                 #3 Number of observations in the data set.
                          as.integer(alpha[i]*n),        #4 Depth of the contour to be calculated.
                          as.double(rep(0,n*(n - 1)/2)), #5 First coordinate of the vertices.
                          as.double(rep(0,n*(n - 1)/2)), #6 Second coordinate of the vertices.
                          as.integer(1),                 #7 Logical signaling empty contour
                          as.integer(1),                 #8 Number of vertices in the contour.
                          as.integer(1),                 #9 Logical signaling whether dithering was performed.
                          PACKAGE = "mrfDepthLight")

      names(Result)[[i]] <- "Contour"
      Result[[i]] <- list(vertices = cbind(TResult[[5]][1:TResult[[8]]],
                                       TResult[[6]][1:TResult[[8]]]
                                       ),
                        depth = TResult[[4]]/n,
                        count = TResult[[8]],
                        empty = as.logical(TResult[[7]]),
                        dithering = as.logical(TResult[[9]]),
                        type = "depthContour")
    }

  Result[[length(alpha) + 1]] <- "hdepth"
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalcHContour <- function(x, alpha, directions, options) {

  NObs <- nrow(x)

  Center <- hdepthmedian(x = x)$median

  Result <- vector("list",length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
    ContResult <- CalcOneHalfContourIntersect(x, Center, alpha[i], directions, options)
    names(Result)[[i]] <- "Contour"
    Result[[i]] <- list(vertices = ContResult$vertices,
                      converged = ContResult$converged,
                      depth = floor(NObs*alpha[i])/NObs,
                      count = nrow(ContResult$vertices),
                      empty = 0,
                      type = "depthContour")
  }

  Result[[length(alpha) + 1]] <- "hdepth"
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalcOneHalfContourIntersect <- function(x, Center, alpha, Directions, options) {

  NObs <- nrow(x)
  NVar <- ncol(x)
  distributionData <- x
  NPoints <- nrow(Directions)

  #Initialise the depth to target
  TargetDepth <- floor(NObs*alpha)/NObs

  #Center the data
  x <- sweep(x, MARGIN = 2, Center, FUN = "-")

  #Find the upper bounds by projection onto the direction.
  #The multivariate depth of points on the contour is smaller
  #than the univariate depth on the projection. Therefore
  #we can take the point with univariate depth equal
  #to the TargetDepth as an upper limit.
  Upperbounds <- matrix(0.0,nrow = NPoints,ncol = NVar)
  for (i in 1:NPoints) {
    ProjectionData <- distributionData %*% Directions[i,]
    #The vector in the matrix Distribution has a direction defined.
    #Therefore the target point must be situated in the second half of the projection.
    TempInd <- NObs - ceiling(TargetDepth*NObs)
    StartDev <- sort(ProjectionData,partial = TempInd)[TempInd]
    Upperbounds[i,] <- 1.2 * abs(StartDev)*Directions[i,]
  }

  #Lets go
  Lowerbounds <- matrix(0.0,ncol = NVar,nrow = NPoints)
  ConsideredPoints <- 0.5 * Upperbounds + Lowerbounds

  #Initialise constant for iterative algorithm
  maxit <- options$maxIter
  count <- 1
  IndToImprove <- rep(TRUE,NPoints)
  DConsPoints <- rep(0.0,NPoints)
  Converged <- rep(TRUE,NPoints)

  #Exception for points equal to the median
  IndCancel <- which(rowSums(Upperbounds) == 0)
  IndToImprove[IndCancel] <- FALSE

  Directions = NULL
  while ((sum(IndToImprove == TRUE) != 0) & (count <= maxit)) {
    Result <- hdepth(x = distributionData,
                     z = ConsideredPoints[IndToImprove,,drop = FALSE],
                     options = options)
    DConsPoints[IndToImprove] <- Result$depthZ
    if(is.null(Directions)){ Directions = Result$directions}
    #Find out which are close enough
    IndStop <- which(abs(DConsPoints - TargetDepth) <= 1/(2*NObs))
    IndToImprove[IndStop] <- FALSE
    DConsPoints[IndStop] <- TargetDepth
    #Find out which need to be closer to the center
    IndCloser <- which(DConsPoints < TargetDepth)
    Upperbounds[IndCloser,] <- ConsideredPoints[IndCloser,]
    ConsideredPoints[IndCloser,] <- 0.5 * (Lowerbounds[IndCloser,] + ConsideredPoints[IndCloser,])
    #Find out which need to be further away from the center
    IndFurther <- which(DConsPoints > TargetDepth)
    Lowerbounds[IndFurther,] <- ConsideredPoints[IndFurther,]
    ConsideredPoints[IndFurther,] <- 0.5 * (ConsideredPoints[IndFurther,] + Upperbounds[IndFurther,])

    #Add a breaking condition for extreme slow convergence situations
    TEMP1 <- sqrt(rowSums(matrix(abs(Lowerbounds - Upperbounds),ncol = NVar)^2))
    TEMP2 <- sqrt(rowSums(Upperbounds^2))
    AidInd <- which(TEMP1/TEMP2 <= 10^-10)
    IndToImprove[AidInd] = FALSE

    count <- count + 1
  }

  ConsideredPoints <- sweep(x = ConsideredPoints, MARGIN = 2, Center, FUN = "+")
  Converged[IndToImprove] <- FALSE

  return(list(vertices = ConsideredPoints,
              converged = Converged)
         )
}

CalcProjTypeContour <- function(x, alpha, type, Directions, options) {

  if (type == "projdepth") {
    TResult <- projdepth(x = x,z = x, options = options)
    Center <- projmedian(x = x, projDepths = TResult$depthZ, options = options)$max
  }
  else if (type == "sprojdepth") {
    TResult <- sprojdepth(x = x,z = x, options = options)
    Center <- sprojmedian(x = x, sprojDepths = TResult$depthZ, options = options)$max
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  Result <- vector("list",length(alpha) + 1)
  names(Result)[[length(alpha) + 1]] <- "type"

  for (i in 1:length(alpha)) {
    ContResult <- CalculateOneProjTypeContour(x, Center, alpha[i], type, Directions, options)
    names(Result)[[i]] <- "Contour"
    Result[[i]] <- list(vertices = ContResult$vertices,
                      converged = ContResult$converged,
                      depth = alpha[i],
                      count = nrow(ContResult$vertices),
                      empty = 0,
                      type = "depthContour")
  }

  Result[[length(alpha) + 1]] <- type
  class(Result) <- c("mrfDepth", "depthContour")
  return(Result)

}

CalculateOneProjTypeContour <- function(x, Center, alpha, type, Directions, options) {

  NDirections <- nrow(Directions)
  NVar <- ncol(x)

  #Center the data
  x <- sweep(x,MARGIN = 2,Center,FUN = "-")

  #Initialise target depth
  TargetDepth <- alpha
  TargetOutlying <- (1/TargetDepth) - 1

  #Find the upper bounds
  Upperbounds <- matrix(0.0,nrow = NDirections, ncol = NVar, byrow = TRUE)
  for (i in 1:NDirections) {
    ProjectionData <- x %*% Directions[i,]
    if (type == "projdepth") {
      Scale <- mad(ProjectionData)
      Location <- median(ProjectionData)
      Upperbounds[i,] <- (TargetOutlying * Scale + Location) * Directions[i,]
    }
    else if (type == "sprojdepth") {
      Quants <- quantile(ProjectionData, probs = c(0.25, 0.5,  0.75), type = 5)
      IQR <- Quants[3] - Quants[1]
      Location <- Quants[2]
      MC <- medcouple(ProjectionData)
      if (MC > 0) { MC <- 3*MC } else{MC <- 4*MC }
      w2 <- Quants[3] + 1.5 * exp(MC) * IQR
      Upperbounds[i,] <- 1.2 * (TargetOutlying * (w2 - Location) + Location )  * Directions[i,]
    }
    else{
      stop("Type not recognized")
    }
  }

  #Find the lower bounds by finding the intersections with the bag.
  Lowerbounds <- matrix(0.0,nrow = NDirections, ncol = NVar, byrow = TRUE)

  # Initialise starting points
  ConsideredPoints <- 0.5 * (Upperbounds + Lowerbounds)

  maxit <- options$maxIter
  count <- 1
  IndToImprove <- rep(TRUE,NDirections)
  DConsPoints <- rep(0.0,NDirections)
  Converged <- rep(TRUE,NDirections)

  while ((sum(IndToImprove == TRUE) != 0) & (count <= maxit)) {
    if (type == "projdepth") {
      DConsPoints[IndToImprove] <- projdepth(x = x,
                                            z = ConsideredPoints[IndToImprove,,drop = FALSE],
                                            options = options)$depthZ
    }
    else if (type == "sprojdepth") {
      DConsPoints[IndToImprove] <- sprojdepth(x = x,
                                             z = ConsideredPoints[IndToImprove,,drop = FALSE],
                                             options = options)$depthZ
    }
    else{
      stop("Type not recognized.")
    }

    #Find out which are close enough
    IndStop <- which(abs(DConsPoints - TargetDepth) <=  10^-3)
    IndToImprove[IndStop] <- FALSE
    DConsPoints[IndStop] <- TargetDepth
    #Find out which need to be closer to the center
    IndCloser <- which(DConsPoints < TargetDepth)
    Upperbounds[IndCloser,] <- ConsideredPoints[IndCloser,]
    ConsideredPoints[IndCloser,] <- 0.5 * (Lowerbounds[IndCloser,] + ConsideredPoints[IndCloser,])
    #Find out which need to be further away from the center
    IndFurther <- which(DConsPoints > TargetDepth)
    Lowerbounds[IndFurther,] <- ConsideredPoints[IndFurther,]
    ConsideredPoints[IndFurther,] <- 0.5 * (ConsideredPoints[IndFurther,] + Upperbounds[IndFurther,])

    #Add a breaking condition for extreme slow convergence situations
    TEMP1 <- sqrt(rowSums(matrix(abs(Lowerbounds - Upperbounds),ncol = NVar)^2))
    TEMP2 <- sqrt(rowSums(Upperbounds^2))
    AidInd <- which(TEMP1/TEMP2 <= 10^-10)
    IndToImprove[AidInd] = FALSE

    count <- count + 1
  }

  ConsideredPoints <- sweep(x = ConsideredPoints, MARGIN = 2, Center, FUN = "+")
  Converged[IndToImprove] <- FALSE

  return(list(vertices = ConsideredPoints,
              converged = Converged)
  )
}
