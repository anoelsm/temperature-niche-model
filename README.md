# Effects of temperature variability and dispersal on realized temperature niches

Computational model in R of phytoplankton competition as a result of overlapping temperature niches across a latitudinal 
gradient of temperature variability and dispersal.

## Description

This repository includes code to run the model (Smith&Barton-model-2023.R), analyze model output (Smith&Barton-model_analysis-2023.R), and plot model results (Smith&Barton-model_figures-2023.R) to recreate figures seen in the above publication.  

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)

## Installation

Project dependencies are listed at the top of the R script as well as a function to install and load the required packages. Additionally, to successfully run the model analysis, users must download and source the following functions: Smith&Barton-model_parameters-2023.R, tniche_analysis-update.R and diversity_function-update.R.

## Usage

Running the code as is, reproduces the results found in the above publication. Additionally, users can explore how changing parameters such as temperature-dependent mortality strength, niche widths, initial conditions, etc. influence the role of dispersal and temperature amplitude on shaping phytoplankton realized niches and diversity. 
