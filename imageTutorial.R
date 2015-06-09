rm(list = ls())

path <- "~/Box Sync/R/iMage001/data/wtAbscess"
setwd(path)

library(devtools)
library(tiff)
library(rtiff)
library(reshape2)

source('~/Box Sync/R/iMage001/iMageFunctions.R')

subdirs <- list.dirs(path, T, F)
subdirs

# This stuff tells us about the images
# (You may need to adjust this for your images)
zstep <- 10        # z-step size in microns
side <- 512        # what is the xy dimension in pixels
pwidth <- 2.4859   # pixel width in microns

# This stuff defines parts of the analysis
# (for now keep these as they are)
npixel <- 5000    # set the number of pixels to be sampled
size <- 30         # set the sampling cube size in pixels

# unique bit of text pertaining to channel 1
ch1 <- "Aa"
# ... and channel 2
ch2 <- "Sg"

# convert the .tif files to an R friendly array
tiffToArray(subdirs = subdirs, 
            side = side, 
            channels = c(ch1, ch2), 
            threshold = F
            )

# where are the new files?
thresholdedStackFiles <- list.files(getwd(), "thresholdedStack", full.names = T)

# do the number crunching on proportion occupancy
proportionOccupancy(files = thresholdedStackFiles, 
                    channels = c(ch1, ch2),
                    zstep = zstep,
                    side = side,
                    pwidth = pwidth,
                    npixel = npixel,
                    size = size
                    )

# read in the results
prod <- read.csv("productivity.csv")
dist <- read.csv("distanceOccupancy.csv")

# plot the productivities
aa_prod <- subset(prod, strain == "Aa")
sg_prod <- subset(prod, strain == "Sg")
plot(aa_prod$proportion ~ sg_prod$proportion,
     ylab = "Aa productivity",
     xlab = "Sg productivity")

# plot the distances
dist$direction <- gsub("ch1", "Aa", dist$direction)
dist$direction <- gsub("ch2", "Sg", dist$direction)

par(mfrow = c(2, 2))
directions <- unique(dist$direction)

for(i in 1:length(directions)){
  s <- subset(dist, direction == directions[i])
  plot(occupancy_mean ~ distance, dist, type = "n", main = directions[i])
  nstacks <- max(s$stack)
  for(i in 1:nstacks){
    s2 <- subset(s, stack == i)
    lines(occupancy_mean ~ distance, s2)
  }
}

