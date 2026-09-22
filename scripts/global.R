#### GLOBAL ####

## packages
library(tidyverse)
library(purrr)
library(ggh4x)
library(cowplot)
library(metR)
library(patchwork)
library(ggimage)
library(png)
library(grid)

## figure themes
theme_set(
  theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size=10),
      legend.position = "bottom",
      legend.title = element_text(size=12),
      legend.text = element_text(size=12),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(size = 10, face = "bold")))

## Ne*s layer
ne_s_layer <- function(){
  geom_function(
    fun = function(x) 1 / (2 * x),
    color = "white",
    linewidth = 0.75,
    linetype = "dashed"
  )
}

## function for saving figures
save_fig<-function(plot, filename, width, height){
  ggsave(
    plot = plot,
    filename = filename,
    device = "tiff",
    width = width*6.5,
    height = height*6.5,
    units = "in",
    scale = 2,
    dpi = 300)
}
