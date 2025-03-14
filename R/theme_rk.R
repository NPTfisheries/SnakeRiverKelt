# Purpose: ggplot2 theme for all figures.
# Author: Ryan N. Kinzer
# Date: 8/25/21

theme_rk <- function(){ 
  font <- "sans"   #assign font family up front
  theme_classic() %+replace%    #replace elements we want to change
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      legend.position = 'bottom',
      #strip.placement = 'outside',
      strip.background = element_rect(colour = NA),
      axis.title=element_text(size=10),
      axis.text = element_text(size=10),
      strip.text = element_text(size=10),
      plot.title = element_text(size = 11, face = 'bold', vjust = 2, hjust = 0),
      plot.subtitle = element_text(size = 10, face = 'plain', vjust = 2, hjust = 0),
      plot.caption = element_text(size = 8, face = 'italic', hjust = 1)
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      # plot.title = element_text(             #title
      #   family = font,            #set font family
      #   size = 20,                #set font size
      #   face = 'bold',            #bold typeface
      #   hjust = 0,                #left align
      #   vjust = 2),               #raise slightly
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      # plot.caption = element_text(           #caption
      #   family = font,            #font family
      #   size = 9,                 #font size
      #   hjust = 1),               #right align
      # axis.title = element_text(             #axis titles
      #   family = font,            #font family
      #   size = 10),               #font size
      # axis.text = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 9),                #font size
      # axis.text.x = element_text(            #margin for axis text
      #   margin=margin(5, b = 10))
    )
}