#!/bin/R
# A helpful R script to plot interspecies comparisons for Wasp paper. 

# Load some libraries:
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
library(cowplot)
library(dplyr)


# Bring in pensylvanica first:  
setwd("~/Desktop/WASPS/vespula_pensylvanica/scaffs/scaffs/")

# Import the scaffolds:
vp.scaff3<-read.table("Scaffold03.pruned.txt", stringsAsFactors = F, header = TRUE)
vp.scaff4<-read.table("Scaffold04.pruned.txt", stringsAsFactors = F, header = TRUE)
vp.scaff5<-read.table("Scaffold05.pruned.txt", stringsAsFactors = F, header = TRUE)
vp.scaff7<-read.table("Scaffold07.pruned.txt", stringsAsFactors = F, header = TRUE)

# Turn off scientific notation
options(scipen = 999)

## SCAFFOLD 03 ## 

vp.scaff3$title <- "Scaffold 3"
vp.Plot3 <- ggplot(data=vp.scaff3, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumpurple4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Plot3

## SCAFFOLD 04 ## 

vp.scaff4$title <- "Scaffold 4"
vp.Plot4 <- ggplot(data=vp.scaff4, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "skyblue4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Plot4

## SCAFFOLD 05 ## 

vp.scaff5$title <- "Scaffold 5"
vp.Plot5 <- ggplot(data=vp.scaff5, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "steelblue2", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Plot5

## SCAFFOLD 07 ## 

vp.scaff7$title <- "Scaffold 7"
vp.Plot7 <- ggplot(data=vp.scaff7, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumseagreen", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Plot7

# Chang working directory to the appropriate directory:
setwd("~/Desktop/WASPS/vespula_pensylvanica/new_recombination/")

# Read in the tripleStats.txt file produced from SlopeCalc.py on HPCC:
vp.recomb3 <- read.table("Scaffold03.recombination.txt", stringsAsFactors = F, header = TRUE)
vp.recomb4 <- read.table("Scaffold04.recombination.txt", stringsAsFactors = F, header = TRUE)
vp.recomb5 <- read.table("Scaffold05.recombination.txt", stringsAsFactors = F, header = TRUE)
vp.recomb7 <- read.table("Scaffold07.recombination.txt", stringsAsFactors = F, header = TRUE)

# Make a list of all the dataframes;
vp.recombination <- list(vp.recomb3, vp.recomb4, vp.recomb5, vp.recomb7)

# Remove NaNs from all the dataframes in the list 
#vp.recombination_clean <- lapply(vp.recombination, function(x) { na.omit(x)})

# Take absolute value of all the slopes: 
#vp.positive_slopes <- lapply(vp.recombination_clean, function(x) {x$absSlope <- abs(x$slope); return(x) })

# Log transform the slope values:
vp.logtransform <- lapply(vp.recombination, function(x) { x$logSlope <- log10(x$slope); return(x) })

# Very small positive values generally turn into larger negative values. 
# Zeros turn into -Inf values. 

names(vp.logtransform) <- c("rcmb3","rcmb4", "rcmb5", "rcmb7")

## RECOMBINATION 03 ## 

vp.logtransform[['rcmb3']]$title <- "Scaffold 3"
vp.Rplot3 <- ggplot(data=vp.logtransform[['rcmb3']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumpurple4")+
  geom_point(colour = "mediumpurple4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Rplot3

## RECOMBINATION 04 ## 

vp.logtransform[['rcmb4']]$title <- "Scaffold 4"
vp.Rplot4 <- ggplot(data=vp.logtransform[['rcmb4']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "skyblue4")+
  geom_point(colour = "skyblue4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Rplot4

## RECOMBINATION 05 ## 

vp.logtransform[['rcmb5']]$title <- "Scaffold 5"
vp.Rplot5 <- ggplot(data=vp.logtransform[['rcmb5']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "steelblue2")+
  geom_point(colour = "steelblue2") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Rplot5

## RECOMBINATION 07 ## 

vp.logtransform[['rcmb7']]$title <- "Scaffold 7"
vp.Rplot7 <- ggplot(data=vp.logtransform[['rcmb7']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumseagreen")+
  geom_point(colour = "mediumseagreen") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vp.Rplot7 

# Now, Consobrina: 

# Set working directory:
setwd("~/Desktop/WASPS/vespula_consobrina/scaffs/")

vc.scaff3 <- read.table("Scaffold03.pruned.txt", stringsAsFactors = F, header = TRUE)
vc.scaff4 <- read.table("Scaffold04.pruned.txt", stringsAsFactors = F, header = TRUE)
vc.scaff5 <- read.table("Scaffold05.pruned.txt", stringsAsFactors = F, header = TRUE)
vc.scaff7 <- read.table("Scaffold07.pruned.txt", stringsAsFactors = F, header = TRUE)

## SCAFFOLD 03 ## 

vc.scaff3$title <- "Scaffold 3"
vc.Plot3 <- ggplot(data=vc.scaff3, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumpurple4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Plot3

## SCAFFOLD 04 ## 

vc.scaff4$title <- "Scaffold 4"
vc.Plot4 <- ggplot(data=vc.scaff4, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "skyblue4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Plot4

## SCAFFOLD 05 ## 

vc.scaff5$title <- "Scaffold 5"
vc.Plot5 <- ggplot(data=vc.scaff5, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "steelblue2", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Plot5

## SCAFFOLD 07 ## 

vc.scaff7$title <- "Scaffold 7"
vc.Plot7 <- ggplot(data=vc.scaff7, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumseagreen", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Plot7

# Set working directory:
setwd("~/Desktop/WASPS/vespula_consobrina/new_recombination/")

# Import the scaffold01 centimorgan x physical distance data:
vc.recomb3 <- read.table("Scaffold03.recombination.txt", stringsAsFactors = F, header = TRUE)
vc.recomb4 <- read.table("Scaffold04.recombination.txt", stringsAsFactors = F, header = TRUE)
vc.recomb5 <- read.table("Scaffold05.recombination.txt", stringsAsFactors = F, header = TRUE)
vc.recomb7 <- read.table("Scaffold07.recombination.txt", stringsAsFactors = F, header = TRUE)

# Make a list of all the dataframes;
vc.recombination <- list(vc.recomb3, vc.recomb4, vc.recomb5, vc.recomb7)

# Remove NaNs from all the dataframes in the list 
#vc.recombination_clean <- lapply(vc.recombination, function(x) { na.omit(x)})

# Take absolute value of all the slopes: 
#vc.positive_slopes <- lapply(vc.recombination_clean, function(x) {x$absSlope <- abs(x$slope); return(x) })

# Log transform the slope values:
vc.logtransform <- lapply(vc.recombination, function(x) { x$logSlope <- log10(x$slope); return(x) })

# Very small positive values generally turn into larger negative values. 
# Zeros turn into -Inf values. 
names(vc.logtransform) <- c("rcmb3","rcmb4", "rcmb5", "rcmb7")

## RECOMBINATION 03 ## 

vc.logtransform[['rcmb3']]$title <- "Scaffold 3"
vc.Rplot3 <- ggplot(data=vc.logtransform[['rcmb3']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumpurple4")+
  geom_point(colour = "mediumpurple4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0), 
                     breaks = c(39e+06, 42e+06, 45e+06, 48e+06, 51e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Rplot3

## RECOMBINATION 04 ## 

vc.logtransform[['rcmb4']]$title <- "Scaffold 4"
vc.Rplot4 <- ggplot(data=vc.logtransform[['rcmb4']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "skyblue4")+
  geom_point(colour = "skyblue4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Rplot4

## RECOMBINATION 05 ## 

vc.logtransform[['rcmb5']]$title <- "Scaffold 5"
vc.Rplot5 <- ggplot(data=vc.logtransform[['rcmb5']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "steelblue2")+
  geom_point(colour = "steelblue2") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Rplot5

## RECOMBINATION 07 ## 

vc.logtransform[['rcmb7']]$title <- "Scaffold 7"
vc.Rplot7 <- ggplot(data=vc.logtransform[['rcmb7']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumseagreen")+
  geom_point(colour = "mediumseagreen") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vc.Rplot7 

# Now Vidua!

# Set working directory:
setwd("~/Desktop/WASPS/vespula_vidua/scaffs/")

# Import the scaffolds:
vv.scaff3 <- read.table("Scaffold03.pruned.txt", stringsAsFactors = F, header = TRUE)
vv.scaff4 <- read.table("Scaffold04.pruned.txt", stringsAsFactors = F, header = TRUE)
vv.scaff5 <- read.table("Scaffold05.pruned.txt", stringsAsFactors = F, header = TRUE)
vv.scaff7 <- read.table("Scaffold07.pruned.txt", stringsAsFactors = F, header = TRUE)

## SCAFFOLD 03 ## 

vv.scaff3$title <- "Scaffold 3"
vv.Plot3 <- ggplot(data=vv.scaff3, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumpurple4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Plot3

## SCAFFOLD 04 ## 

vv.scaff4$title <- "Scaffold 4"
vv.Plot4 <- ggplot(data=vv.scaff4, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "skyblue4", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Plot4

## SCAFFOLD 05 ## 

vv.scaff5$title <- "Scaffold 5"
vv.Plot5 <- ggplot(data=vv.scaff5, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "steelblue2", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Plot5

## SCAFFOLD 07 ## 

vv.scaff7$title <- "Scaffold 7"
vv.Plot7 <- ggplot(data=vv.scaff7, aes(x=Position, y=centiMorgan, group=1)) +
  geom_point(colour = "mediumseagreen", size=1) + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6, accuracy = 1.0),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Plot7

# Chang working directory to the appropriate directory:
setwd("~/Desktop/WASPS/vespula_vidua/new_recombination/melee_calc/")

# Read in the tripleStats.txt file produced from SlopeCalc.py on HPCC:
vv.recomb3 <- read.table("Scaffold03.recombination.txt", stringsAsFactors = F, header = TRUE)
vv.recomb4 <- read.table("Scaffold04.recombination.txt", stringsAsFactors = F, header = TRUE)
vv.recomb5 <- read.table("Scaffold05.recombination.txt", stringsAsFactors = F, header = TRUE)
vv.recomb7 <- read.table("Scaffold07.recombination.txt", stringsAsFactors = F, header = TRUE)

# Make a list of all the dataframes;
vv.recombination <- list(vv.recomb3, vv.recomb4, vv.recomb5, vv.recomb7)

# Remove NaNs from all the dataframes in the list 
#vv.recombination_clean <- lapply(vv.recombination, function(x) { na.omit(x)})

# Take absolute value of all the slopes: 
#vv.positive_slopes <- lapply(vv.recombination_clean, function(x) {x$absSlope <- abs(x$slope); return(x) })

# Log transform the slope values:
vv.logtransform <- lapply(vv.recombination, function(x) { x$logSlope <- log10(x$slope); return(x) })

# Very small positive values generally turn into larger negative values. 
# Zeros turn into -Inf values. 
names(vv.logtransform) <- c("rcmb3","rcmb4", "rcmb5", "rcmb7")
## RECOMBINATION 03 ## 

vv.logtransform[['rcmb3']]$title <- "Scaffold 3"
vv.Rplot3 <- ggplot(data=vv.logtransform[['rcmb3']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumpurple4")+
  geom_point(colour = "mediumpurple4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6),
                     breaks = c(39e+06, 42e+06, 45e+06, 48e+06, 51e+06), limits = c(38.5e+06,51.1e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Rplot3

## RECOMBINATION 04 ## 

vv.logtransform[['rcmb4']]$title <- "Scaffold 4"
vv.Rplot4 <- ggplot(data=vv.logtransform[['rcmb4']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "skyblue4")+
  geom_point(colour = "skyblue4") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6),
                     breaks = c(52e+06, 54e+06, 56e+06, 58e+06, 60e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Rplot4

## RECOMBINATION 05 ## 

vv.logtransform[['rcmb5']]$title <- "Scaffold 5"
vv.Rplot5 <- ggplot(data=vv.logtransform[['rcmb5']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "steelblue2")+
  geom_point(colour = "steelblue2") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6),
                     breaks = c(62e+06, 64e+06, 66e+06, 68e+06, 70e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Rplot5

## RECOMBINATION 07 ## 

vv.logtransform[['rcmb7']]$title <- "Scaffold 7"
vv.Rplot7 <- ggplot(data=vv.logtransform[['rcmb7']], aes(x=midpoint, y=logSlope, group=1)) +
  geom_line(colour = "mediumseagreen")+
  geom_point(colour = "mediumseagreen") + theme_bw() +
  scale_x_continuous(labels = scales::unit_format(unit = NULL, scale = 1e-6),
                     breaks = c(79e+06, 81e+06, 83e+06, 85e+06, 87e+06)) +
  theme(axis.title = element_text(face="plain",size=15)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + facet_grid(. ~ title)
vv.Rplot7 

# For 24 line plots, in a 4x6 grid:
grid <- plot_grid(vp.Plot3, vp.Plot4, vp.Plot5, vp.Plot7, 
          vp.Rplot3, vp.Rplot4, vp.Rplot5, vp.Rplot7,
          vc.Plot3, vc.Plot4, vc.Plot5, vc.Plot7, 
          vc.Rplot3, vc.Rplot4, vc.Rplot5, vc.Rplot7,
          vv.Plot3, vv.Plot4, vv.Plot5, vv.Plot7, 
          vv.Rplot3, vv.Rplot4, vv.Rplot5, vv.Rplot7, 
          ncol = 4, nrow = 6, labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
                                       "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"),
          align = "v") 

grid.arrange(grid, left = textGrob("Recombination rate (-log(cM/Mb))", rot = 90),
             bottom = textGrob("physical position (Mb)"))


# Set working directory:
setwd("~/Desktop/WASPS/final_figures/")

# Set the high-resolution png file 
png(filename = "~/Desktop/WASPS/Vespula.comparison.recombination.map.png", 
    width = 15, height = 15, res = 400, units = 'in', 
    type = c("quartz")) # for high resolution figure 
dev.off()


