
# Load packages.
library(geodata)
library(sf)
library(tidyverse)
library(terra)

# Make an object to help navigate the subdirectories.
my_dir_path <- getwd()
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$figs    <- file.path( my_dir_path, "figs" )
wd$data    <- file.path( my_dir_path, "data" )

# Check for presence of subdirectories. Create if needed.
invisible({
  lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) )
})

proj_eqd    <- "+proj=aeqd +lon_0=-100 +lat_0=40 +datum=WGS84 +units=m +no_defs"
proj4_wgs   <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj4_tpers <- "+proj=tpers +h=15000000 +azi=5 +lon_0=-90 +lat_0=40"
