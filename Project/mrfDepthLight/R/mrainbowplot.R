mrainbowplot <- function(x, depths, col=NULL, LegendTitle="") {

  ######
  # Check input.
  if (missing(x)) { stop("Input argument x is required.") }
  if (missing(depths)) { stop("Input argument depths is required.") }

  #Check x
  x <- data.matrix(x)
  if (!is.matrix(x)) { stop("x must be a matrix.") }
  if ( sum(is.nan(x)) != 0 ) { stop("x contains missing cases.") }
  n1 <- dim(x)[1]
  p1 <- dim(x)[2]
  if (p1 != 2) { stop("Only implemented for bivariate data.") }
  PlotData <- data.frame(cbind(x,depths))
  colnames(PlotData) <- c("x","y","depth")
  if (is.null(colnames(x))) {
    Labs <- c("variable 1", "variable 2")
  } else {
    Labs <- colnames(x)
  }


  #Check depths
  depths <- data.matrix(depths)
  if ( sum(is.nan(depths)) != 0 ) { stop("depths contains missing cases.") }
  n2 <- dim(depths)[1]
  p2 <- dim(depths)[2]
  if (n1 != n2) { stop("depths does not have the same size as x.")}
  if (p2 != 1) { stop("depths must be column matrix.")}

  #Check rgb colors
  if (is.null(col)) { col <- makeColors_MRainbow() }
  col <- data.matrix(col)
  if ( sum(is.nan(col)) != 0 ) { stop("col contains missing cases.") }
  p3 <- dim(col)[2]
  if (p3 != 3) { stop("col must have three columns.")}
  if (sum(col > 1) != 0) { stop("All values in the paramter col must lie in [0,1].") }
  if (sum(col < 0) != 0) { stop("All values in the paramter col must lie in [0,1].") }
  RGBCols <- rgb(col[,1],col[,2],col[,3])


  #Create basic plot
  plot <- ggplot()
  plot <- plot + geom_point(data = PlotData,
                          mapping = aes_string(x = "x", y = "y",
                                                  colour = "depth"),
                          size = 4
                             )
  #Make a colorbal scale
  plot <- plot + scale_colour_gradientn(colours = RGBCols, name = LegendTitle)
  plot

  #Set up the figure
  plot <- plot + ggtitle( "Bivariate rainbowplot. " )
  xRange <- extendrange(x[,1], f = 0.05)
  yRange <- extendrange(x[,2], f = 0.05)
  plot <- plot + coord_cartesian(xlim = xRange, ylim = yRange)

  #give plot the package look
  plot <- plot + mrfDepth_theme() +
    scale_x_continuous(name = Labs[1]) +
    scale_y_continuous(name = Labs[2])

  return(plot)

}

makeColors_MRainbow <- function() {
  RGBmatrix <- c(0.823529411764706, 1.000000000000000, 0.823529411764706,
                0.392156862745098, 0.823529411764706, 0.392156862745098,
                0.235294117647059, 0.784313725490196, 0.235294117647059,
                0.117647058823529, 0.666666666666667, 0.117647058823529,
                0.058823529411765, 0.509803921568627, 0.058823529411765,
                0                , 0.196078431372549, 0)
  RGBmatrix <- matrix(RGBmatrix,ncol = 3,byrow = TRUE)
}
