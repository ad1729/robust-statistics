bagdistance <- function(x, z = NULL, options = list()){

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) { stop("The input argument x must be a numeric data matrix.") }
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
  # check options
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
  if ( "approx"  %in% names(options) ) {
    approx <- options[["approx"]]
    if (!is.logical(approx)) {
      stop("Option approx must be a logical value.")
    }
  } else {
    options$approx <- TRUE
  }

  #####
  # Check data for possible exact fit situations.
  tol <- 1e-7
  ScaledX <- scale(x)
  temp <- attributes(ScaledX)
  ColumnSd <- temp[["scaled:scale"]]
  if (sum(ColumnSd <= 1e-14) > 0) {
    warning("One of the variables has zero standard deviation zero. Check the data matrix x.")
    return(list(bagdistance = NULL,
                flag = NULL,
                converged = NULL,
                dimension = sum(ColumnSd > 1e-14),
                hyperplane = as.numeric(ColumnSd <= 1e-14),
                type = "bagdistance")
           )
  }
  w1 <- try(svd(ScaledX/sqrt(n1 - 1)),silent = TRUE)
  if (!is.list(w1)) {
    warning("The singular-value decomposition of the data matrix x could not be computed.")
    return(list(bagdistance = NULL,
                flag = NULL,
                converged = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = "bagdistance")
           )
  }
  if (min(w1$d) < tol) {
    warning("An exact fit is found. Check the output for more details.")
    return(list(bagdistance = NULL,
                flag = NULL,
                converged = NULL,
                dimension = sum(w1$d > tol),
                hyperplane = w1$v[,which(w1$d == min(w1$d))[1]],
                type = "bagdistance")
           )
  }
  #####
  # Prepare the start of the algorithm
  #####
  # Find the median of the depths of the data points.
  Result1 <- hdepth(x = x, z = x, options = options)
  if (is.null(Result1$depthZ)) {
    stop("The halfspace depth of x can not be computed.")
  }
  # Find the halfspace median of x.
  Result2 <- hdepthmedian(x = x)
  if (is.null(Result2$median)) {
    stop("The halfspace median of x can not be computed.")
  }
  Center <- Result2$median

  # Initialise constants.
  NVar <- p1
  NObs <- n1
  # Note that in the hdepth call z was taken to be x.
  TargetDepth <- sort(Result1$depthZ, decreasing = TRUE)[ceiling(n1/2)]
  NPoints <- n2

  # Initialise vector of halfspace distances.
  Distances <- rep(0.0,NPoints)

  #####
  # Start the real calculations.
  #####
  if (NVar == 1) {
    # This case is pretty easy.
    Depths <- Result1$depthX
    Ind <- which(Depths >= TargetDepth)
    LowerBound <- min(x[Ind])
    UpperBound <- max(x[Ind])

    Ind <- which(z <= Center)
    Distances[Ind] <- abs(z[Ind] - Center) / abs(LowerBound - Center)

    Ind <- which(z >= Center)
    Distances[Ind] <- abs(z[Ind] - Center) / abs(UpperBound - Center)
    cutoff <- sqrt(qchisq(0.99,p1))
    FlagZ <- Distances <= cutoff

    return(list(bagdistance = Distances,
                flag = FlagZ,
                converged = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = "bagdistance"))

  }
  else if ( (NVar == 2) & (options$approx==FALSE) ) {
    # First calculate the isodepth contour and then calculate the intersection of each ray with the contour.

    # Center all data
    distributionData <- sweep(x,MARGIN = 2, Center, FUN = "-")
    calcDistData <- sweep(z,MARGIN = 2, Center, FUN = "-")

    # Compute the isodepth contour and calculate intersections.
    VerticesContour <- depthContour(distributionData,
                                    alpha = TargetDepth,
                                    type = "hdepth"
                                    )$Contour
    if ("dataDimension" %in% names(VerticesContour)) {
      # An exact fit was found
      stop("The isodepth contour could not be computed.")
    }
    VerticesContour <- VerticesContour$vertices
    NVertices <- nrow(VerticesContour)

    # Define an angle transformation function such that all points have increasing angles in [0, 2pi)
    AngleTransfo <- function(Angles,Data){
      Ind180 <- which(Data[,1] < 0)
      Ind360 <- which(Data[,2] < 0 & Data[,1] > 0)

      Angles[Ind180] <- Angles[Ind180] + pi
      Angles[Ind360] <- Angles[Ind360] + 2*pi
      return(Angles)
    }

    # Find the angels of all points.
    ContAngles <- atan(VerticesContour[,2]/VerticesContour[,1])
    ContAngles <- AngleTransfo(Angles = ContAngles, Data = VerticesContour)

    DistriAngles <- atan(distributionData[,2]/distributionData[,1])
    DistriAngles <- AngleTransfo(Angles = DistriAngles, Data = distributionData)

    CalcDistAngles <- atan(calcDistData[,2]/calcDistData[,1])
    CalcDistAngles <- AngleTransfo(Angles = CalcDistAngles, Data = calcDistData)

    # Sort all angles keeping track of which group the angles belong to.
    TAngles <- rbind(cbind(ContAngles,rep(1.0,length(ContAngles))),
                    cbind(CalcDistAngles,rep(0.0,NPoints))
                    )
    TData <- rbind(VerticesContour,calcDistData)
    TInd <- order(TAngles[,1],TAngles[,2],decreasing = FALSE)
    TAngles <- TAngles[TInd,]
    TData <- TData[TInd,]

    # Start calculating the intersections.
    DistriInd <- which(TAngles[,2] == 1)
    NDistriInd <- length(DistriInd)

    # Do calculations for points between 2pi and 0.
    ToCalcInd <- c()
    if (DistriInd[1] != 1) {ToCalcInd <- 1:(DistriInd[1] - 1)}
    if (DistriInd[NDistriInd] != length(TInd)) {ToCalcInd <- c(ToCalcInd,(DistriInd[NDistriInd] + 1):length(TInd))}
    if (!is.null(ToCalcInd)) {
      # Find the equation of the line element in the contour.
      if (TData[DistriInd[1],1] == TData[DistriInd[NDistriInd],1]) {
        # The line segment is vertical
        Temp <- TData[ToCalcInd,]
        if (is.vector(Temp)) {Temp <- matrix(Temp,ncol = 2)}
        Intersections <- cbind(rep(TData[DistriInd[1],1], (DistriInd[1] - 1)),
                              (Temp[,2]/Temp[,1])*TData[DistriInd[1],1]
        )
        Distances[TInd[ToCalcInd]-NVertices] <- sqrt(Temp[,1]^2 + Temp[,2]^2) / sqrt(Intersections[,1]^2 + Intersections[,2]^2)
      }
      else{
        Temp <- TData[ToCalcInd,]
        if (is.vector(Temp)) {Temp <- matrix(Temp,ncol = 2)}
        SegmentSlope <- (TData[DistriInd[1],2] - TData[DistriInd[NDistriInd],2])/(TData[DistriInd[1],1] - TData[DistriInd[NDistriInd],1])
        IntersectX <- (SegmentSlope * TData[DistriInd[1],1] - TData[DistriInd[1],2]) / (SegmentSlope - Temp[,2]/Temp[,1])
        IntersectY <- (Temp[,2]/Temp[,1])*IntersectX
        Intersections <- cbind(IntersectX,IntersectY)
        Distances[TInd[ToCalcInd]-NVertices] <- sqrt(Temp[,1]^2 + Temp[,2]^2) / sqrt(Intersections[,1]^2 + Intersections[,2]^2)
      }
    }
    for (i in 1:(length(DistriInd) - 1)) {
      if ((DistriInd[i] + 1) != DistriInd[i + 1]) {
        ToCalcInd <- (DistriInd[i] + 1):(DistriInd[i + 1] - 1)
        # Find the equation of the line element in the contour.
        if (TData[DistriInd[i],1] == TData[DistriInd[i + 1],1]) {
          # The line segment is vertical
          Temp <- TData[ToCalcInd,]
          if (is.vector(Temp)) {Temp <- matrix(Temp,ncol = 2)}
          Intersections <- cbind(rep(TData[DistriInd[i],1],nrow(Temp)),
                                (Temp[,2]/Temp[,1])*TData[DistriInd[i],1]
          )
          Distances[TInd[ToCalcInd] - NVertices] <- sqrt(Temp[,1]^2 + Temp[,2]^2) / sqrt(Intersections[,1]^2 + Intersections[,2]^2)
        }
        else{
          Temp <- TData[ToCalcInd,]
          if (is.vector(Temp)) {Temp <- matrix(Temp,ncol = 2)}
          SegmentSlope <- (TData[DistriInd[i],2] - TData[DistriInd[i + 1],2])/(TData[DistriInd[i],1] - TData[DistriInd[i + 1],1])
          IntersectX <- (SegmentSlope * TData[DistriInd[i],1] - TData[DistriInd[i],2]) / (SegmentSlope - Temp[,2]/Temp[,1])
          IntersectY <- (Temp[,2]/Temp[,1])*IntersectX
          Intersections <- cbind(IntersectX,IntersectY)
          DistTemp <- sqrt(Temp[,1]^2 + Temp[,2]^2) / sqrt(Intersections[,1]^2 + Intersections[,2]^2)
          Distances[TInd[ToCalcInd] - NVertices] <- DistTemp
        }
      }
    }
    Distances[is.nan(Distances)] <- 0

    cutoff <- sqrt(qchisq(0.99,p1))
    FlagZ <- Distances <= cutoff

    return(list(bagdistance = Distances,
                flag = FlagZ,
                converged = NULL,
                dimension = NULL,
                hyperplane = NULL,
                type = "bagdistance"))
  }
  else{

    # Center all data
    distributionData <- sweep(x,MARGIN = 2, Center, FUN = "-")
    calcDistData <- sweep(z,MARGIN = 2, Center, FUN = "-")

    # Get the directions
    Directions <- calcDistData
    Directions <-  Directions/ sqrt(rowSums(Directions^2))

    Result <-  CalcOneHalfContourIntersect(x = distributionData,
                                              Center = rep(0.0,NVar),
                                              alpha = TargetDepth,
                                              Directions = Directions,
                                              options = options
                                              )
    Intersects <- Result$vertices

    Distances <- sqrt( rowSums(calcDistData^2) )
    Scales <- sqrt( rowSums(Intersects^2) )
    ScaledDistances <- Distances / Scales
    ScaledDistances[is.nan(ScaledDistances)] <- 0

    cutoff <- sqrt(qchisq(0.99,p1))
    FlagZ <- ScaledDistances <= cutoff

    return(list(bagdistance = ScaledDistances,
                flag = FlagZ,
                converged = Result$converged,
                dimension = NULL,
                hyperplane = NULL,
                type = "bagdistance"
                )
    )

    }
}
