sprojdepth <- function(x, z = NULL, options = NULL) {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  #Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n1 <- nrow(x)
  p1 <- ncol(x)
  if (n1 > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  #Check the z data.
  if (is.null(z)) { z <- x }
  z <- data.matrix(z)
  if (!is.numeric(z)) { stop("The input argument z must be a numeric data matrix.") }
  n2 <- nrow(z)
  p2 <- ncol(z)
  if (p1 != p2) { stop("Data dimension has to be the same for x and z.") }
  if (n2 > sum(complete.cases(z))) { stop("Missing values in z are not allowed.") }
  #check options
  if (is.null(options)) { options <- list() }
  if (!is.list(options)) { stop("options must be a list") }

  #####
  #Check data for possible exact fit situations.
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
  w1 <- try(svd(ScaledX/sqrt(n1 - 1)),silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x could not be computed.")
    return(list(depthX = NULL,
                depthZ = NULL,
                singularSubsets = NULL,
                cutoff = NULL,
                flagX  =  NULL,
                flagZ = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = "spdepth")
    )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit was found. Check output for more details.")
    return(list(depthX = NULL,
                depthZ = NULL,
                singularSubsets  =  NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[which(w1$d == min(w1$d))[1]],
                type = "spdepth")
    )
  }

  Original <- options(warn = 1)
  Result <- adjOutlyingness(x = x, z = z, options = options)
  options(warn = Original$warn)

  if (!is.null(Result$hyperplane)) {
    return(list(depthX = NULL,
                depthZ = NULL,
                singularSubsets = NULL,
                cutoff = NULL,
                flagX = NULL,
                flagZ = NULL,
                dimension = NULL,
                hyperplane = Result$hyperplane,
                SubspaceIndicator = Result$SubspaceIndicator,
                type = "sprojdepth"
                )
           )
  } else{
    return(list(depthX = 1/(1 + Result$outlyingnessX),
                depthZ = 1/(1 + Result$outlyingnessZ),
                singularSubsets = Result$singularsubsets,
                cutoff = 1/(1 + Result$cutoff),
                flagX = Result$flagX,
                flagZ = Result$flagZ,
                type = "sprojdepth")
           )
  }

}
