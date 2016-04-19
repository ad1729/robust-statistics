bagplot <- function(CompBagResult,
                    colorbag = NULL,
                    colorloop = NULL,
                    colorchull = NULL,
                    databag = TRUE,
                    dataloop = TRUE,
                    plotFence = FALSE){
  ######
  # Check input.
  if (missing(CompBagResult)) { stop("Input argument CompBagResult is required.") }

  if ( !("mrfDepth" %in% class(CompBagResult)) ) { stop("The input parameter CompBagResult should be the return of a call to compBagplot. ") }
  if ( !("compBagplot" %in% class(CompBagResult)) ) { stop("The input parameter CompBagResult should be the return of a call to compBagplot. ") }
  if (is.null(colorbag)) { colorbag = rgb(0.6,0.6,1) }
  if (is.null(colorloop)) { colorloop = rgb(0.8,0.8,1) }
  if (is.null(colorchull)) { colorchull = rgb(1,1,1) }

  IndBag <- which(CompBagResult$datatype[,3] == 1)
  dataContBag <- data.frame(CompBagResult$bag)
  dataBag <- data.frame(CompBagResult$datatype[IndBag,1:2,drop = FALSE])
  IndLoop <- which(CompBagResult$datatype[,3] == 2)
  dataLoop <- data.frame(CompBagResult$datatype[IndLoop,1:2,drop = FALSE])
  dataFence <- data.frame(CompBagResult$fence)
  IndOutl <- which(CompBagResult$datatype[,3] == 3)
  dataOutl <- data.frame(CompBagResult$datatype[IndOutl,1:2,drop = FALSE])
  colnames(dataFence) <- c("x","y")
  colnames(dataLoop) <- c("x","y")
  colnames(dataBag) <- c("x","y")
  colnames(dataContBag) <- c("x","y")
  colnames(dataOutl) <- c("x","y")

  XLabel <- colnames(CompBagResult$datatype)[1]
  YLabel <- colnames(CompBagResult$datatype)[2]

  # Initialise plot
  plot <- ggplot()

  # Add the fence
  if (plotFence) {
    dataReg <- dataFence[chull(dataFence),]
    colnames(dataReg) <- c("x","y")
    dataReg <- rbind(dataReg, dataReg[1,])
    plot <- plot + geom_path(data = dataReg,
                             mapping = aes_string(x = "x", y = "y"),linetype = "dashed"
                             )
  }

  # Add the filled loop
  dataReg <- dataLoop[chull(dataLoop),]
  colnames(dataReg) <- c("x","y")
  plot <- plot + geom_polygon(data = dataReg,
                             mapping = aes_string(x = "x", y = "y"),
                             fill = colorloop)

  # Add the filled bag
  dataReg <- dataContBag[chull(dataContBag),]
  colnames(dataReg) <- c("x","y")
  plot <- plot + geom_polygon(data = dataReg,
                             mapping = aes_string(x = "x", y = "y"),
                             fill = colorbag)

  # Add the innermost contour if applicable
  if (nrow(CompBagResult$chull) > 1) {
    dataCHull <- data.frame(CompBagResult$chull)
    colnames(dataCHull) <- c("x","y")
    plot <- plot + geom_polygon(data = dataCHull,
                               mapping = aes_string(x = "x", y = "y"),
                               fill = colorchull)
  }

  # Add the center
  dataReg <- data.frame(matrix(CompBagResult$center,ncol = 2))
  colnames(dataReg) <- c("x","y")
  plot <- plot + geom_point(data = dataReg,
                           mapping = aes_string(x = "x", y = "y"),
                           shape = 23, size = 4, color = "red", fill = "red"
                           )
  # If requested plot the points
  if (databag) {
    plot <- plot + geom_point(data = dataBag,
                             mapping = aes_string(x = "x", y = "y")
                    )
  }
  if (dataloop) {
    plot <- plot + geom_point(data = dataLoop,
                             mapping = aes_string(x = "x", y = "y")
                             )
  }

  # Plot outliers
  if (nrow(dataOutl) > 0) {
    plot <- plot + geom_point(data = dataOutl,
                             mapping = aes_string(x = "x", y = "y"),
                             shape = 8,color = "red"
                             )
  }

  # Set up the figure
  PlotData <- CompBagResult$datatype[,-3]
  colnames(PlotData) <- c("x","y")
  PlotData <- rbind(PlotData, dataFence)
  xRange <- extendrange(PlotData[,1], f = 0.05)
  yRange <- extendrange(PlotData[,2], f = 0.05)
  plot <- plot + coord_cartesian(xlim = xRange, ylim = yRange)

  # give plot the package look
  plot <- plot + mrfDepth_theme()

  # Finalise
  plot <- plot + xlab(XLabel) + ylab(YLabel)
  if (CompBagResult$type == "hdepth") {
    TitleLabel <- paste("Bagplot based on halfspace depth")
  } else if (CompBagResult$type == "projdepth") {
    TitleLabel <- paste("Bagplot based on projection depth")
  } else if (CompBagResult$type == "sprojdepth") {
    TitleLabel <- paste("Bagplot based on skew-adjusted projection depth")
  } else{
    TitleLabel <- paste("Bagplot based on", CompBagResult$type, "depth")
  }
  plot <- plot + ggtitle( TitleLabel )

  return(plot)

}
