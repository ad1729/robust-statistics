sdepth <- function(x, z=NULL) {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  #Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  if (p1 > 2) { stop("Simplicial depth can only be computed for p<=2.") }
  if (n1 < (p1 + 1)) {stop("At least (p+1) observations are needed.") }

  #Check z argument
  if (is.null(z)) { z <- x }
  z <- data.matrix(z)
  if (!is.numeric(z)) { stop("The input argument z must be a numeric data matrix.") }
  n2 <- nrow(z)
  p2 <- ncol(z)
  if (p1 != p2) { stop("The data dimension has to be the same for x and z.") }
  if (n2 > sum(complete.cases(z))) { stop("Missing values in z are not allowed.") }

  #####
  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  ScaledX <- scale(x)
  temp <- attributes(ScaledX)
  ColumnSd <- temp[["scaled:scale"]]
  if (sum(ColumnSd <= 1e-14) > 0) {
    warning("One of the variables has zero standard deviation. Check the data matrix x.")
    return(list(depthZ = NULL,
                dimension = sum(ColumnSd > 1e-14),
                hyperplane = as.numeric(ColumnSd <= 1e-14),
                type = "sdepth"
    )
    )
  }
  w1 <- try(svd(ScaledX/sqrt(n1 - 1)),silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x could not be computed.")
    return(list(depthZ = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = "sdepth"
    )
    )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check the output for more details.")
    return(list(depthZ = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
                type = "sdepth"
    )
    )
  }


  if (p1 == 1) {
    depth <- rep(NaN,n2)
    TotLines <- (n1*(n1 - 1))/2
    for (i in 1:n2) {
      NSmaller <- sum(x <= z[i])
      NBigger <- sum(x >= z[i])
      depth[i] <- (NSmaller*NBigger)/TotLines
    }
    ReturnedResult <- list(depthZ = depth,
                           dimension = NULL,
                           hyperplane = NULL,
                           type = "sdepth")
  }
  else if (p1 == 2) {
    Result <- .Fortran("HSDEP2",
                       as.double(z[,1,drop = TRUE]), #1 First coordinate of points to compute the depth of.
                       as.double(z[,2,drop = TRUE]), #2 Second coordinate of points to compute the depth of.
                       as.integer(n2),               #3 Number of points in the set z.
                       as.double(x[,1,drop = TRUE]), #4 First coordinate of the data matrix x.
                       as.double(x[,2,drop = TRUE]), #5 Second coordinate of the data matrix x.
                       as.integer(n1),               #6 Number of reference points.
                       as.double(rep(0,n2)),         #7 Halfspace depth of the points in z.
                       as.double(rep(0,n2)),         #8 Simplicial depth of the points in z.
                       PACKAGE = "mrfDepthLight")

    ReturnedResult <- list(depthZ = Result[[8]],
                           dimension = NULL,
                           hyperplane = NULL,
                           type = "sdepth")

  }
  return(ReturnedResult)
}
