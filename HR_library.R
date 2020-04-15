# Description ####
# Functions for analysing and constructing developmental trees

# Written by B. Spanjaard, 2018

# Dependencies ####
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(collapsibleTree))
suppressPackageStartupMessages(require(data.tree))

# Parameters ####
theme_bs <- 
  theme_update(text = element_text(size = 24),
               panel.background = element_rect(fill = "white", colour = "black"),
               panel.grid.major.y = element_line(colour = "grey"),
               panel.grid.major.x = element_line(colour = "grey"),
               legend.key = element_rect(fill = "white"),
               plot.title = element_text(hjust = 0.5))

# Functions ####
