# Load required libraries
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)

# load data
envData <- USE.MCMC::Worldclim_tmp
envData <- terra::rast(envData, type="xyz")
