adjOutlyingness <- function(x, z=NULL, options = list() ){

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {	stop("The input argument x must be a numeric data matrix.") }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  # Check the z data.
  if (is.null(z)) { z <- x }
  z <- data.matrix(z)
  if (!is.numeric(z)) { stop("The input argument z must be a numeric data matrix.") }
  n2 <- nrow(z)
  p2 <- ncol(z)
  if (p1 != p2) { stop("The data dimension has to be the same for x and z.") }
  if (n2 > sum(complete.cases(z))) { stop("Missing values in z are not allowed.") }

  #check options en load defaults if necessary.
  if (is.null(options)) { options = list()  }
  if (!is.list(options)) { stop("options must be a list") }
  if ( "type" %in% names(options) ) {
    type <- options[["type"]]
  } else {
    type <- "Affine"
  }
  if ( "ndir"  %in% names(options) ) {
    ndir <- options[["ndir"]]
  } else {
    ndir <- NULL
  }
  if ( "seed"  %in% names(options) ) {
    seed <- options[["seed"]]
  } else {
    seed <- NULL
  }

  # Check number of directions and type.
  typeID <- match(type,c("Affine","Rotation","Shift"))[1]
  if (is.na(typeID)) { stop("The input parameter type must be one of: Affine, Rotation or Shift.") }
  if (is.null(ndir)) {
    if (typeID == 1)	ndir <- 250*p1
    if (typeID == 2)	ndir <- 250*20
    if (typeID == 3)	ndir <- 250*50
  }
  CalcAll <- 0
  # If the specified number of directions is larger than the
  # possible different directions, switch to exact computation.
  if (is.numeric(ndir)) {
    if (ndir < 1) { stop("The number of directions must be a positive integer.") }
    if (typeID == 1) {
      ndir0 <- choose(n1,p1)
      if (ndir0 <= ndir) {
        ndir <- ndir0
        CalcAll <- 1
      }
    }
    if (typeID == 2) {
      ndir0 <- choose(n1,2)
      if (ndir0 <= ndir) {
        ndir <- ndir0
        CalcAll <- 1
      }
    }
  }
  if (!is.numeric(ndir)) {
    if (ndir == "all") {
      if (typeID == 1) {
        ndir <- choose(n1,p1)
        CalcAll <- 1
        if (ndir > 1e7) { stop("ndir is larger than 1e7. Try a smaller value of ndir.") }
      }
      if (typeID == 2) {
        ndir <- choose(n1,2)
        CalcAll <- 1
        if ( ndir > 1e7 ) { stop("ndir is larger than 1e7. Try a smaller value of ndir.") }
      }
      if (typeID == 3) { stop("Cannot compute all directions for type Shift.") }
    }
    else stop("The input parameter ndir is not recognized.")
  }
  if ( (n1 < (p1 + 1)) & typeID == 1) { stop("When type is affine, n should be larger than p.") }
  if (is.null(seed)) {
    seed <- 10
  }
  if (!is.numeric(seed)) {
    stop("The seed must be a strictly positive integer.")
  }
  if (seed <= 0) {
    stop("The seed must be a strictly positive integer.")
  }



  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  ScaledX <- scale(x)
  temp <- attributes(ScaledX)
  ColumnSd <- temp[["scaled:scale"]]
  if (sum(ColumnSd <= 1e-14) > 0) {
    warning("One of the variables has zero standard deviation. Check the data matrix x.")
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = sum(ColumnSd > 1e-14),
                hyperplane = as.numeric(ColumnSd <= 1e-14),
                inSubspace = NULL,
                type = "adjOutlyingness")
           )
  }
  w1 <- try(svd(ScaledX/sqrt(n1 - 1)),silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x could not be computed.")
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = NULL,
                hyperplane = NULL,
                inSubspace = NULL,
                type = "adjOutlyingness")
    )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit is found. Check the output for more details.")
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
                inSubspace = NULL,
                type = "adjOutlyingness")
    )
  }

  #####
  # Perform the actual computations
  x <- rbind(x,z)
  n <- nrow(x)
	p <- ncol(x)

	Result <- .C("adjprojout",
		         as.integer(n),		    	#1 Total number of points.
		         as.integer(p),		    	#2 Dimension of the data.
  		       as.integer(ndir),      #3 Number of directions.
  		       as.double(x),		    	#4 Data matrix (both x and z).
	  	       as.double(rep(0,n)),  	#5 Computed adjusted outlyingness.
		         as.double(0),	   		  #6 Medcouple computed on the adj outl. values of x.
		         as.integer(0),		  	  #7 Number of singular directions.
		         as.integer(typeID), 	  #8 Integer indicating which type of directions to consider.
		         as.integer(n1),			  #9 Number of observations in x.
		         as.integer(CalcAll),   #10 Flag indicating whether all possible directions should be considered.
		         as.double(rep(0,p1)),  #11 Vector containing the direction on which a zero IQR is found.
		         as.integer(seed),      #12 The seed.
		         PACKAGE = "mrfDepthLight")

	AdjOutlyingness <- Result[[5]]
  if (sum(abs(Result[[11]])) > tol) {
    warning("A direction was found for which the robust scale estimate zero. See the help page for more details.", call. = FALSE)
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = NULL,
                hyperplane = Result[[11]],
                inSubspace = as.logical(AdjOutlyingness),
                type = "adjOutlyingness")
    )
  }

  cutoff <- sqrt(qchisq(0.99,p1)) * median(AdjOutlyingness[1:n1])

  FlagX <- AdjOutlyingness[1:n1] <= cutoff
  FlagZ <- AdjOutlyingness[(n1 + 1):(n1 + n2)] <= cutoff

  ReturnedResult <- list(outlyingnessX = AdjOutlyingness[1:n1],
                         outlyingnessZ = AdjOutlyingness[(n1 + 1):(n1 + n2)],
                         cutoff = cutoff,
                         flagX = FlagX,
                         flagZ = FlagZ,
                         type = "adjOutlyingness",
                         singularsubsets = Result[[7]],
                         dimension = NULL,
                         hyperplane = NULL,
                         inSubspace = NULL,
                         type = "adjOutlyingness")
	class(ReturnedResult) <- "mrfDepth"
	return(ReturnedResult)
}
