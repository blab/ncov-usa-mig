#File: color_schemes.R
#Author(s): Amin Bemanian
#Date: 07/07/24
#Description: Color scheme definitions for all the figures in paper, figures should generally use these colors,
# and any future Shiny or d3 applications should also be based on the same colorways
#Arguments: None

library(RColorBrewer)

#Blue to Orange gradient for RR heatmaps and chloropleths, log based scale
#Defaults of 1.1 and 2.0
RR_log_grad <- function(LB,UB){
    scale_fill_gradient2(name="RR",
                        high = "#D67C34",
                        low = "#4C90C0",
                        limits=c(log10(LB),log10(UB)),
                        breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                        labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB))
}