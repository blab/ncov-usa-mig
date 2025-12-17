#File: color_schemes.R
#Author(s): Amin Bemanian
#Date: 07/07/24
#Description: Color scheme definitions for all the figures in paper, figures should generally use these colors,
# and any future Shiny or d3 applications should also be based on the same colorways
#Arguments: None

library(RColorBrewer)

##Region scales
REGION_SCALE <- c("Far West" = "#faa", 
                  "Great Lakes" = "#e66", 
                  "Mideast" = "#ea5",
                  "New England"="#fdb",
                  "Plains"="#fd0", 
                  "Rocky Mountain"="#fea",
                  "Southeast"="#ad4",
                  "Southwest"="#bfb", 
                  "Mexico"="#0fd",
                  "Eastern Canada"="#9df", 
                  "Western Canada"="#aaf")
region_fill_scale <- function(){
  scale_fill_manual(values=REGION_SCALE,
                    name="BEA Regions",
                    na.translate = FALSE,
                    na.value = "grey80")
}

region_color_scale <- function(){
  scale_color_manual(values=REGION_SCALE,
                     name="BEA Regions",
                     na.translate = FALSE,
                     na.value = "grey80")
}

##Country scales
COUNTRY_SCALE <- c("Mexico" = "#5d4",
                   "USA" = "#35d",
                   "Canada" = "#d46")
country_fill_scale <- function(){
  scale_fill_manual(values=COUNTRY_SCALE,
                    name="Country",
                    na.translate = FALSE,
                    na.value = "grey80")
}

country_color_scale <- function(){
  scale_color_manual(values=COUNTRY_SCALE,
                     name="Country",
                     na.translate = FALSE,
                     na.value = "grey80")
}


RR_HIGH <- "#D67C34"
RR_LOW <- "#4C90C0"
#Blue to Orange gradient for RR heatmaps and chloropleths, log based scale
#Defaults of 1.1 and 2.0
RR_log_grad <- function(LB=0.5,UB=2.0){
    scale_fill_gradient2(name="RR",
                        high = RR_HIGH,
                        low = RR_LOW,
                        limits=c(log10(LB),log10(UB)),
                        breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                        labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB))
}

fill_bound <- function(x,LB=0.5,UB=2.0){
  log_x <- log10(x)
  pmax(pmin(log_x, log10(UB)), log10(LB))
}


#Manuallly feed breaks
RR_log_grad_manual <- function(b){
  l = length(b)
  scale_fill_gradient2(name="RR",
                       high = RR_HIGH,
                       low = RR_LOW,
                       limits=c(log10(min(b)),
                                log10(max(b))),
                       breaks =c(log10(b)),
                       labels = c(b))
}

#Bounds should be provided in linear terms
log_seq_scale <- function(LB=1E3,UB=1E6){
  if(LB >= UB){
    stop("Upper bound must be larger than lower bound!")
  }
  color_vec = c("#eee","#fe6","#e51")
  log_LB = log10(LB)
  log_UB = log10(UB)
  log_diff = log_UB - log_LB
  
  
  log_breaks = round(c(log_LB,
                 log_LB + log_diff * 1/3,
                 log_LB + log_diff * 2/3,
                 log_UB))
  label_expr <- parse(text = paste0("10^", log_breaks))
  print(label_expr)
  
  scale_fill_gradientn(name="Number of Sequences",
                       colors = color_vec,
                       limits=c(log_LB,log_UB),
                       breaks =log_breaks,
                       labels = label_expr)
}

#Bounds should be provided in linear terms
effort_scale <- function(LB=0,UB=5000,step_size=1000){
  if(LB >= UB){
    stop("Upper bound must be larger than lower bound!")
  }
  color_vec <- c("#eee","#fe6","#e51")
  eff_breaks <- seq(from=LB, to=UB, by=step_size)
  eff_labels <- formatC(eff_breaks, format = "d", big.mark = ",")
  eff_labels[length(eff_labels)] <- paste0(eff_labels[length(eff_labels)], "+")
  
  
  scale_fill_gradientn(name="Sequences per 100,000",
                       colors = color_vec,
                       limits=c(LB,UB),
                       breaks = eff_breaks,
                       labels = eff_labels)
}


## Sex categorical scales
SEX_COLORS <- c(
  "Female" = "#395",
  "Male" = "#C4C"
)
sex_fill_scale <- function(){
  scale_fill_manual(values=SEX_COLORS,na.translate=FALSE, name = "Sex")
}
sex_color_scale <- function(){
  scale_color_manual(values=SEX_COLORS,na.translate=FALSE, name = "Sex")
}

