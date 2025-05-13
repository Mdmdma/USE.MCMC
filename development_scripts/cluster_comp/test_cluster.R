# Load required libraries
Sys.setenv(https_proxy="http://proxy.service.consul:3128")
Sys.setenv(http_proxy="http://proxy.service.consul:3128")
library(USE.MCMC)
library(terra)
library(virtualspecies)
library(sf)
library(ggplot2)
library(parallel)

cat("loading of libraries works")
