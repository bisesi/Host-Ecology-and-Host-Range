# Bacteriophage diet breadth is impacted by interactions between bacteria

![Bacterial ecology](https://github.com/bisesi/Host-Ecology-and-Host-Range/figures/final-figs/imgs/figure-1.png)

### Authors:

Ave T. Bisesi, Wolfram Moebius, Carey Nadell, Eleanore G. Hansen, Steven D. Bowden, William R. Harcombe

## Repository Overview

This code repository contains the necessary data and scripts to reproduce the figures and analyses presented in the manuscript, "Bacteriophage diet breadth is impacted by interactions between bacteria." 

To reproduce the analyses, clone or download the project and follow the instructions provided below.

## Getting started

1. Clone or download the repository to your local machine.
2. If you are using R Studio, open the R project file, ```Host-Ecology-and-Host-Range.Rproj```. This will set up your working directories correctly, so there is no need to change them manually. 

## Figure generation

Figures can be recreated using the individual scripts found in the folder ```figures/final-figs/scripts```. The script ```generate-final-imgs.R``` will generate the data for each figure and export the figure to the correct location as a png. Statistics for Figures 4, 5 and Supplemental Figure 3 can be generated and saved as xlsx files using the ```generate-stats-tables.R``` found in the same folder.  

## Experimental data

Tecan and PFU data for all relevant experiments can be found in the ```experimental-data/tecan-data``` folder. Each experiment is labelled by date; the README in each folder provides details about the conditions tested. Helper functions to clean and analyze the data can be found in the ```functions``` folder. 

## Required R Packages

Before running the scripts, ensure you have the following R packages installed:

- ```tidyverse```: A collection of R packages designed for data science. 
- ```ggpubr```: Functions for creating publication-ready plots.
- ```rstatix```: Functions for performing various statistical tests.
- ```deSolve```: Functions for the numerical treatment of systems of differential equations.
- ```ODESensitivity```: Functions for sensitivity analysis in ordinary differential equation (ode) models.
- ```cowplot```: Functions for creating publication-ready plots.
- ```ggtext```: Functions for creating publication-ready plots.
- ```zoo```: Functions for methods for totally ordered indexed observations.
- ```patchwork```: Functions for creating publication-ready plots.
- ```readxl```: Functions to load xlsx files.