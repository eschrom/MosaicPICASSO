# Mosaic-PICASSO
Mosaic-PICASSO, an R Package

## Overview
Mosaic-PICASSO is an algorithm that infers and removes artefacts due to autofluorescence and spillover from multiplex images of biological tissues. The PICASSO algorithm was developed by Seo et al. (https://doi.org/10.1038/s41467-022-30168-z) and the Mosaic addition by Cang et al. (https://doi.org/10.1093/bioinformatics/btad784). This R package implementation includes a few extra features, including the option to use correlation instead of mutual information as an internal metric, the option to split channels into mutually exclusive groups for independent correction, and the option to enhance contrast after correction. 

## Contents
In addition to the source R code, this package also includes:
- documentation detailing the purpose, inputs, parameters, and outputs of each function
- a tutorial in vignette and html formats describing how to use Mosaic-PICASSO

## Installation
The Mosaic-PICASSO R package can be installed from Github:
- Install the R package "devtools" if you don't already have it.
- Run the command devtools::install_github("eschrom/MosaicPICASSO", build_vignettes = TRUE)
- Restart your R session

Installation should only take a few minutes, but perhaps more if many of the dependencies must also be installed or updated.

## Computing Requirements and Runtimes
This implementation of MosaicPICASSO was developed on a PC laptop with a 11th Gen Intel(R) Core(TM) i9-11900H (2.50 GHz) processor and 16.0 GB RAM, running Windows 10 64-bit operating system, and using R version 4.3.1 "Beagle Scouts." Computing times scale with the size of the image and the number of channels; in general, computing times will be minutes to hours.
