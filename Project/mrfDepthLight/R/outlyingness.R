outlyingness <- function(x,z=NULL, options=list()) {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
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
  if (is.null(options)) { options <- list()  }
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
  if ( "stand"  %in% names(options) ) {
    stand <- options[["stand"]]
  } else {
    stand <- "MedMad"
  }
  if ( "centered"  %in% names(options) ) {
    centered <- options[["centered"]]
  } else {
    centered <- FALSE
  }
  if ( "h"  %in% names(options) ) {
    h <- options[["h"]]
  } else {
    h <- NULL
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
    if (typeID == 1) ndir <- 250*p1
    if (typeID == 2) ndir <- 250*20
    if (typeID == 3) ndir <- 250*50
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
        if (ndir > 1e7) { stop("ndir is larger than 1e7. Try a smaller value of ndir.") }
      }
      if (typeID == 3) { stop("Cannot compute all directions for type Shift.") }
    } else stop("The input parameter ndir is not recognized.")
  }
  if (p1 == 1)	ndir <- 1
  if (n1 < (p1 + 1) & typeID == 1) { stop("When type is affine, n should be larger than p.") }
  # Check stand argument
  scaleID <- match(stand,c("MedMad","unimcd"))[1]
  if (is.na(scaleID)) { stop("The input parameter stand should be either MedMad or unimcd.") }
  # Check centered option
  if ((centered %in% c(FALSE,TRUE)) == FALSE) {  stop("The input parameter centered should be TRUE or FALSE.") }
  # C code expects the opposite of documented option
  centered <- !centered
  # Check input h.
  if (is.null(h)) { h <- floor(n1/2) + 1 }
  if (!is.numeric(h)) { stop("The input parameter h should be numeric.") }
  if (h < (floor(n1/2) + 1) | h > n1 )	{ stop("The input parameter h should lie between [(n/2)]+1 and n.") }
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
  w1 <- try(svd(scale(x)/sqrt(n1 - 1)),silent = TRUE)
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
                type = "outlyingness")
    )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
                inSubspace = NULL,
                type = "outlyingness")
    )
  }

  #####
  #Do the actual computations
  x <- rbind(x,z)
  n <- nrow(x)
  p <- ncol(x)

	Factor <- ifelse(h<n1,qchisq(h/n1,df = 1),1)

  Result <- .C("projoutlyingness",
          			as.integer(n),	    	 #1 Total number of points.
          			as.integer(p1),	    	 #2 Dimension of the data.
          			as.integer(ndir),	     #3 Number of directions.
          			as.double(x),		       #4 Data matrix (both x and z).
          			as.double(rep(0,n)),	 #5 Computed Stahel-Donoho Outlyingness.
          			as.integer(0),	       #6 Number of singular directions.
          			as.integer(typeID),	   #7 Integer indicating which type of directions to consider.
          			as.integer(n1),	    	 #8 Number of points in reference set x.
          			as.integer(scaleID),	 #9 Integer indicating which type of directions to consider.
          			as.integer(h),	       #10 Integer indicating how many points should be included in the h-subset of unimcd.
          			as.integer(centered),  #11 Integer indicating whether data should be centered.
          			as.double(Factor),     #12 The reweighting factor for unimcd.
          			as.integer(CalcAll),   #13 Flag indicating whether all possible directions should be considered.
          			as.double(rep(0,p1)),  #14 Vector containing the direction on which a zero spread is found.
          			as.integer(seed),      #15 The seed.
          			PACKAGE = "mrfDepthLight")

  Outlyingness <- Result[[5]]
  cutoff <- sqrt(qchisq(0.99,p1)) * median(Outlyingness[1:n1])
  FlagX <- (Outlyingness[1:n1] <= cutoff)
  FlagZ <- (Outlyingness[(n1 + 1):(n1 + n2)] <= cutoff)

  if (sum(abs(Result[[14]])) > tol) {
    warning("A direction was found for which the robust scale equals zero. See the help page for more details.", call. = FALSE)
    return(list(outlyingnessX = NULL,
                outlyingnessZ = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                singularSubsets = NULL,
                dimension = NULL,
                hyperplane = Result[[14]],
                inSubspace = as.logical(Outlyingness),
                type = "outlyingness")
    )
  }

  ReturnedResult <- list(outlyingnessX = Outlyingness[1:n1],
                         outlyingnessZ = Outlyingness[(n1 + 1):(n1 + n2)],
                         cutoff = cutoff,
                         flagX = FlagX,
                         flagZ = FlagZ,
                         singularSubsets = Result[[6]],
                         dimension = NULL,
                         hyperplane = NULL,
                         inSubspace = NULL,
                         type = "outlyingness")

  class(ReturnedResult) <- "mrfDepth"
  return(ReturnedResult)
}
