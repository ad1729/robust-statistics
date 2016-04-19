compBagplot <- function(x, type="hdepth", sizesubset=500, extraDirections = FALSE, options=NULL){

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n <- nrow(x)
  p <- ncol(x)
  if (n > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  if (p != 2) { stop("The bagplot can only be drawn for bivariate data.") }
  if (n < 10) { stop("At least 10 data points are required.") }
  if (is.null(colnames(x))) { colnames(x) <- c("variable 1", "variable 2") }

  # Check type
  Indtype <- match(type,c("hdepth","projdepth","sprojdepth"))[1]
  if (is.na(Indtype)) { stop("The input parameter type should be one of hdepth, projdepth or sprojdepth.") }

  # Check sizesubset
  if (!is.numeric(sizesubset)) { stop("The input parameter sizesubset should be numeric.") }

  # Check ExtraDirections
  if (!is.logical(extraDirections)) { stop("The input parameter ExtraDirections should be a logical.") }


  #check options
  if (is.null(options)) { options <- list()  }
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



  #####
  # Check data for possible exact fit situations.
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
    return(list(center = NULL,
                chull = NULL,
                bag = NULL,
                datatype = NULL,
                flag = NULL,
                depth = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = type)
           )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    return(list(center = NULL,
                chull = NULL,
                bag = NULL,
                datatype = NULL,
                flag = NULL,
                depth = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
                type = type)
           )
  }

  if (type == "hdepth") {
    result <- CompHalfSpaceBagplot(x,sizesubset)
  }
  else if (type == "projdepth" || type == "sprojdepth") {
    result <- CompProjectionTypeBagplot(x, type, extraDirections, options)
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  return(result)

}

CompHalfSpaceBagplot <- function(x,sizesubset) {
  n <- nrow(x)
  p <- ncol(x)
  interpol <- matrix(0.0,2*n,2)
  datatyp <- matrix(0.0,n,3)
  datatyp2 <- matrix(0.0,n,2)
  chull <- matrix(0.0,sizesubset*(sizesubset - 1)/2,2)
  pxpy <- matrix(0.0,n,3)
  PDepths <- matrix(0.0,n,2)
  storage.mode(interpol) <- "double"
  storage.mode(datatyp) <- "double"
  storage.mode(datatyp2) <- "double"
  storage.mode(chull) <- "double"
  storage.mode(pxpy) <- "double"
  storage.mode(PDepths) <- "double"

  result <- .Fortran("bagplot",
     			          as.integer(n),			 	 #1	Number of points
     			          as.double(x[,1]),		 	 #2	First variable of x
     			          as.double(x[,2]),			 #3	Second variable of x
               			as.integer(3),		     #4	Specify whisker format
               			as.integer(sizesubset),#5	nsub
               			as.double(rep(0,2)),	 #6	tukm
               			as.integer(1),				 #7	Number of points in chull
               			chull,					       #8	Contour points of deepest region
               			interpol,				       #9	Points in the bag
               			as.integer(0),				 #10	Number of points in bag
               			datatyp,				       #11	Data with their type
               			as.integer(rep(0,n)),  #12	Indices of outliers
               			as.integer(0),				 #13	Number of outliers
               			datatyp2,				       #14	datatyp2
               			pxpy,					         #15	Aid variable
               			as.integer(1),				 #16	Flag indicating if half of the points ly on a vertical line
               			as.integer(1),				 #17	Flag to indicate whether dithering was performed
               			PDepths,				       #18 Depths of points
               			PACKAGE = "mrfDepthLight")

	if (result[[13]] > 0) {
   	IndOutl <- result[[12]][1:result[[13]]]
  } else {
    IndOutl <- NULL
  }
  Flag <- rep(0,n)
  Flag[IndOutl] <- 1

  if (result[[16]] == 1) {
   		stop("At least half of the observations lie on a line, a univariate boxplot should be used.")
  }

  Datatyp <- result[[11]]
  colnames(Datatyp) <- c(colnames(x), "positionIndicator")

  #Calculate the fence.
  Bag <- result[[9]][1:result[[10]],]
  Center <- result[[6]]
  CenteredBag <- sweep(Bag,MARGIN = 2,Center,FUN = "-")
  Fence <-  3*CenteredBag
  Fence <- sweep(Fence,MARGIN = 2,Center,FUN = "+")

  Result <- list(center = result[[6]],
                chull = result[[8]][1:result[[7]],],
                bag = result[[9]][1:result[[10]],],
                fence = Fence,
                datatype = Datatyp,
                flag = Flag,
                depth = result[[18]][,1],
                type = "hdepth"
                )
  class(Result) <- c("mrfDepth", "compBagplot")

  return(Result)
}

CompProjectionTypeBagplot <- function(x, type, extraDirections, options){

  if (type == "projdepth") {
    Result <- projdepth(x = x,z = x, options = options)
    Center <- projmedian(x = x, projDepths = Result$depthZ)$max
  }
  else if (type == "sprojdepth") {
    Result <- sprojdepth(x = x,z = x, options = options)
    Center <- sprojmedian(x = x, sprojDepths = Result$depthZ)$max
  }
  else{
    stop("The input parameter type is not recognized.")
  }

  n <- nrow(x)
  p <- ncol(x)
  nBag <- ceiling(n/2)
  IndOutl <- (Result$depthZ < Result$cutoff)

  SDepth <- sort(Result$depthZ,decreasing = TRUE)
  BagCutoff <- SDepth[nBag]
  IndBag <- which(Result$depthZ >= BagCutoff)
  Data1 <- x[IndBag,1:2]
  TInd <- chull(Data1[,1],Data1[,2])
  Bag <- Data1[TInd,]

  Datatyp <- cbind(x, rep(0.0,n))
  Datatyp[!IndOutl,3] <- 2
  Datatyp[IndBag,3] <-  1
  Datatyp[IndOutl,3] <- 3

  colnames(Datatyp) <- c(colnames(x), "positionIndicator")

  TInd <- which(Datatyp[,3] != 3)
  IndLoop <- chull(x = Datatyp[TInd,1], y = Datatyp[TInd,2])
  Loop <- Datatyp[TInd[IndLoop],1:2 ]
  Directions <- sweep(rbind(Loop, Bag),MARGIN = 2,Center,FUN = "-")
  if (extraDirections) {
    Angle <- seq(0,2*pi,length.out = 250)
    Directions <- rbind(Directions, cbind(cos(Angle),sin(Angle)))
  }
  Directions <- Directions / sqrt(rowSums(Directions^2))

  Fence <- CalculateOneProjTypeContour(x = x,
                                      Center = Center,
                                      alpha = Result$cutoff,
                                      type = type,
                                      Directions = Directions,
                                      options = options)
  Fence <- Fence$vertices

  result <- list(center = matrix(Center,ncol = p),
                 chull = matrix(Center,ncol = p),
                 bag = Bag,
                 fence = Fence,
                 datatype = Datatyp,
                 flag = as.numeric(Result$flagZ),
                 depth = Result$depthZ,
                 type = type)

  class(result) <- c("mrfDepth", "compBagplot")
  return(result)

}

