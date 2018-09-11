#' FUNCTIONS TO MAP AREA ESTIMATES AND THEIR UNCERTAINTY

#' Function to generate the individual plots of areas, confidence intervals and
#' margins of error
#' 
#' @param totalarea: Numerical value indicating the total valid area
#' @param xlabels: Vector of n labels to be applied to the x axis (e.g. years)
#' @param area: Vector of n area values for a given class
#' @param lower: Vector of n lower confidence interval area values for a given class
#' @param upper: Vector of n upper confidence interval area values for a given class
#' @param mappedarea: Vector of n mapped area values for a given class
#' @param me: Vector of n margin of error values for a given class
#' @param miny: Numerical value indicating min y axis value
#' @param maxy: Numerical value indicating max y axis value
#' @param plot_title: String with plot title
#' @param plotmode: Numerical value indicating type of plot to produce. Values are:
#' 1 for plot with no connecting lines for periods where confidence interval crosses 0,
#' 2 for dotted connected lines for periods wheere the confidence interval crosses 0 and
#' 3 for similar plots to option 1 but with larger labels, suitable for individually
#' saved plots. 
#' @return List with area plot and margin of error plot for a given class.

plot_areas = function(totalarea, xlabels, area, lower, upper, mappedarea, me, miny, maxy, plot_title, plotmode){
  # Need two copies bc of the complexity of the graph
  tempdf = as.data.frame(cbind(seq(1,length(xlabels)), area, lower, upper, mappedarea, me))
  names(tempdf) = c("Years", "Area", "Lower", "Upper", "Mapped_area", "Margin_error")
  tempdf2 = tempdf
  
  # Try to format with padded zeros if xlabels are single digits
  xlabels = sprintf("%02i",xlabels)
  
  # Find rows where the CI or the area go below 0 and assign NA's
  ind_area = which(tempdf$Area < 0)
  ind_lower = which(tempdf$Lower < 0)
  
  tempdf[union(ind_area, ind_lower), 'Area'] = NA
  tempdf[ind_lower, 'Lower'] = NA
  tempdf[ind_lower, 'Upper'] = NA
  
  # Define breaks manually so that the two y axes grids match
  bks = seq(miny, maxy, length.out = 6)
  bks2 = bks / totalarea * 100
  
  # Global options for all plots, unless overrriden below
  yscale = scale_y_continuous(labels=function(n){format(n, scientific = FALSE, big.mark = ",")},
                              breaks=bks, limits=c(miny,maxy), 
                              sec.axis = sec_axis(~./totalarea * 100, breaks=bks2, 
                                                  labels=function(n){format(n, digits=2)},
                                                  name="Percentage of total area\n"))
  ribbon = geom_ribbon(data=tempdf, aes(x=Years, ymin=Lower, ymax=Upper), fill="deepskyblue4", alpha=0.3, size=0.5)
  markers = geom_point(data=tempdf, aes(x=Years, y=Area), shape=3, size=1, stroke=0.5)
  map_area = geom_line(data=tempdf2, aes(x=Years, y=Mapped_area), colour="red", size=0.2)
  
  # Remove CI and area when it intersects with zero
  if (plotmode == 1){
    
    lowerline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Lower), linetype=8, size=0.2)
    upperline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Upper), linetype=8, size=0.2) 
    centerline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Area), linetype=8, size=0.2)
    
    # Or keep the original outlines but remove the fill
  } else if (plotmode == 2) {
    
    lowerline = geom_line(data=tempdf2, aes(x=Years, y=Lower), linetype=8)
    upperline = geom_line(data=tempdf2, aes(x=Years, y=Upper), linetype=8) 
    centerline = geom_line(data=tempdf2, aes(x=Years, y=Area), linetype=8)
    
    # Or keep plotmode1 format but make labels and markers bigger
  } else if (plotmode == 3) {
    
    lowerline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Lower), linetype=8, size=0.5)
    upperline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Upper), linetype=8, size=0.5) 
    centerline = geom_line(data=tempdf[!is.na(tempdf$Area),], aes(x=Years, y=Area), linetype=8, size=0.5)
    markers = geom_point(data=tempdf, aes(x=Years, y=Area), shape=3, size=1, stroke=0.5)
    map_area = geom_line(data=tempdf2, aes(x=Years, y=Mapped_area), colour="red", size=0.5) 
  }
  
  
  # Plot areas with CI. Add dashed lines on top of ribbon for places where CI < 0
  regular_theme = theme(plot.title = element_text(size=7), axis.text=element_text(size=7),
                        axis.title=element_text(size=7)) 
  
  area_plot = ggplot() +  
    lowerline + upperline + centerline + ribbon + map_area + markers +
    scale_x_continuous(breaks=seq(1,length(xlabels)), labels=xlabels, minor_breaks = NULL) +
    yscale + ylab("Area and 95% CI [kha]\n") + xlab("\nTime") +
    ggtitle(plot_title)  + geom_hline(yintercept = 0, size=0.1) + regular_theme
  
  # Plot margin of error
  me_plot = ggplot(data=tempdf, aes(x=Years, y=Margin_error * 100)) + geom_line(size=1.1) + 
    scale_x_continuous(breaks=seq(1,length(xlabels)), labels=xlabels, minor_breaks = NULL) + ylim(0, 200) + 
    ylab("Margin of error [%]\n") + xlab("\nTime") + ggtitle(plot_title) + regular_theme
  
  # Change fontsize if labels biglabels is False
  small_theme = theme(axis.title=element_blank(), 
                      text = element_text(size=7, family="Times New Roman"),
                      plot.title = element_text(size=7, family="Times New Roman"),
                      axis.text.x=element_text(size=7, family="Times New Roman"), 
                      axis.text.y=element_text(size=7, family="Times New Roman"),
                      axis.ticks = element_line(size = 0.1),
                      panel.grid = element_line(size=0.2)) 
  area_plot_small = area_plot + small_theme
  me_plot_small = me_plot + small_theme
  
  if (plotmode == 1 | plotmode == 2){
    return_list = list(area_plot_small, me_plot_small)
  } else if (plotmode == 3) {
    return_list = list(area_plot, me_plot)
  }
  
  return(return_list)
}
