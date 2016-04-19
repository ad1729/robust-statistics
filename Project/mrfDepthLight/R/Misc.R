
mrfDepth_theme <- function(){
  A <- theme_classic(base_size = 12, base_family = "") +
    theme(panel.background =  element_rect(fill = NA, colour = "black", size = 0.25),
          plot.margin = ggplot2::unit(c(1,1,1,1), "cm")
    )
  return(A)
}

GridPlot <- function(plotlist=NULL, VarLayout=NULL) {
  numPlots <- length(plotlist)
  NRow <- nrow(VarLayout)
  NCol <- ncol(VarLayout)

  if (numPlots == 1) {
    print(plotlist[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(NRow, NCol)))

    # Make each plot, in the correct location
    for (i in 1:NRow) {
      for (j in 1:NCol) {
        if (VarLayout[i,j] != 0)
        print(plotlist[[ VarLayout[i,j] ]] , vp = viewport(layout.pos.row = i,
                                                           layout.pos.col = j))
      }
    }
  }
}

matSubstract.c <- compiler:::cmpfun(function(mat, Center, NRow){
  mat - rep(1,NRow) %*% Center
})
