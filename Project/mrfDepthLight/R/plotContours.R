plotContours <- function(x, depthContour, data = TRUE) {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }
  if (missing(depthContour)) { stop("Input argument depthContour is required.") }

  # Check the x data.
  x <- data.matrix(x)
  if (!is.numeric(x)) {  stop("The input argument x must be a numeric data matrix.") }
  n <- nrow(x)
  p <- ncol(x)
  if (n > sum(complete.cases(x))) { stop("Missing values in x are not allowed.") }
  if (p != 2) { stop("The data matrix x must be two dimensional.") }
  if (is.null(colnames(x))) { colnames(x) <- c("variable 1", "variable 2") }

  if ( !("mrfDepth" %in% class(depthContour)) ) { stop("depthContour should be the retun of a call to depthContour.") }
  if ( !("depthContour" %in% class(depthContour)) ) { stop("depthContour should be the retun of a call to depthContour.") }

  # Initialise plot
  plot <- ggplot()

  # If requested plot the points
  xData <- data.frame(x)
  colnames(xData) <- c("x", "y")
  if (data) {
    plot <- plot + geom_point(data = xData,
                              mapping = aes_string(x = "x", y = "y"),
                              col = "grey"
                              )
  }

  #Add the contours
  for (i in 1:(length(depthContour) - 1) ) {
    TResult <- depthContour[[i]]
    if ( TResult$empty != 1) {
      # Add the filled bag
      Ind <- chull(TResult$vertices)
      if(length(Ind)>1){
        dataCont <- data.frame(TResult$vertices[chull(TResult$vertices),])
        dataCont <- rbind(dataCont, dataCont[1,])
        colnames(dataCont) <- c("x","y")
        plot <- plot + geom_path(data = dataCont,
                                 mapping = aes_string(x = "x", y = "y")
        )
      }
    }
  }


  # Set up the figure
  # xRange <- extendrange(x[,1], f = 0.05)
  # yRange <- extendrange(x[,2], f = 0.05)
  # plot <- plot + coord_cartesian(xlim = xRange, ylim = yRange)

  # give plot the package look
  plot <- plot + mrfDepth_theme()

  # Finalise
  plot <- plot + xlab(colnames(x)[1]) + ylab(colnames(x)[2])
  if (depthContour$type == "hdepth") {
    TitleLabel <- paste("Halfspace depth contours")
  } else if (depthContour$type == "projdepth") {
    TitleLabel <- paste("Projection depth contours")
  } else if (depthContour$type == "sprojdepth") {
    TitleLabel <- paste("Skew-adjusted projection depth contours")
  } else{
    TitleLabel <- ""
  }
  plot <- plot + ggtitle( TitleLabel )


  return(plot)

}
