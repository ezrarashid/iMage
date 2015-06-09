### functions for iMage
# image analysis for microbial communities

# 1. thresholdImages()
# which converts the images to binary and stores them in the parent directory as arrays
# if threshold = T then this function also performs a threshold
tiffToArray <- function(
  subdirs, 
  side, 
  channels = NULL,
  threshold = T
) {
  for(f in 1:length(subdirs)){
    # list the file paths
    files <- dir(subdirs[f], ".tif", full.names = T)
    if(is.null(channels) == F){
      # split into channels
      cFiles <- list()
    }
    for(j in 1:length(channels)){
      cFiles[[j]] <- files[grep(channels[j], files)]
      # load em up
      cArray <- array(0, c(side, side, length(cFiles[[j]])))
      for(i in 1:length(cFiles[[j]])){
        print.noquote(paste("loading ch", j, "  ", i, "/", 
                      length(cFiles[[j]]), cFiles[[j]][i]))
        cArray[,,i] <- readTIFF(cFiles[[j]][i])
      }
      # loop to threshold a stack
      c_t <- array(0, c(side, side, length(cFiles[[j]])))
      for(i in 1:dim(cArray)[3]){
        if(threshold == T){
          print.noquote(paste("thresholding ch", j, "  ", i, "/", 
                        length(cFiles[[j]]), cFiles[[j]][i]))
          if(sum(cArray[,,i]) > 0){
            th <- autoThreshold(cArray[,,i], mean(cArray[,,i]))[2]
          }else{th <- 0}
          c_t[,,i][cArray[,,i] > th] <- 1  
        }else
          if(threshold == F){
            c_t <- cArray
          }
      }
      saveRDS(c_t, file = paste(subdirs[f], "_", channels[j], "_thresholdedStack.R", sep = ""))
    }
  }
# end function  
}


#--------------------------------------------------------------------------------------------------


# 2. proportionOccupancy()
# which calculates for n pixels and given distance windows what proportion of the avaiable space withing a distance window is occupied by a given strain.

proportionOccupancy <- function(
  files,
  channels,
  zstep,
  side,
  pwidth,
  size = 30,
  npixel = 5000,
  focus = c("ch1", "ch2"),
  target = c("ch1", "ch2")
) {
  
  print.noquote("Setting up...")
  
  stopifnot(size%%zstep == 0)
  
  # reset npixels
  npixels <- npixel
  
  #### 1. GENERATE A NULL BOX
  
  # generate a null box
  null_box <- expand.grid(x = seq((-size*pwidth), (size*pwidth), by = pwidth), 
                          y = seq((-size*pwidth), (size*pwidth), by = pwidth), 
                          z = seq(-size, size, by = zstep))
  # calculate the distances from focal pixel to each pixel in box
  d <- 0
  center <- null_box[null_box$x == 0 & null_box$y == 0 & null_box$z == 0,]
  for(j in 1:dim(null_box)[1]){
    d[j] <- sqrt((center$x - null_box$x[j])^2 +
                 (center$y - null_box$y[j])^2 +
                 (center$z - null_box$z[j])^2)
  }
  # bins of distances
  ds <- seq(0, (sqrt((size^2)*3))*pwidth, by = 2)
  # for each bin create a vector of positions
  positions <- list()
  for(i in 1:length(ds)){
    positions[[i]] <- which(d <= ds[i] & d > ds[i]-1)
  }

  # objects to store results in
  productivities <- data.frame(file = NULL, channel = NULL, productivity = NULL)
  distance_results <- data.frame(file = NULL, distance = NULL, direction = NULL, 
                                 propOcc_mean = NULL, propOcc_sd = NULL)
  
  #### 2. LOAD THE THRESHOLDED IMAGES
  ch1_files <- files[grep(channels[1], files)]
  ch2_files <- files[grep(channels[2], files)]
  
  for(k in 1:length(ch1_files)){
    
    print.noquote(paste("stack number: ", k, "/", length(ch1_files)))
    print.noquote("    loading...")
    # load em up
    ch1_t <- readRDS(ch1_files[k])
    ch2_t <- readRDS(ch2_files[k])

    # proportion occupied
    r1 <- length(which(ch1_t > 0))/length(ch1_t)
    r2 <- length(which(ch2_t > 0))/length(ch2_t)
    
    theseProductivities <- cbind(k, c(channels[1], channels[2]), c(r1, r2))
    productivities <- rbind(productivities, theseProductivities)
    
    #### 3. GENERATE DISTANCE/OCCUPANCY DATA
    
    # Calculations
    # 1. randomly sample a pixel
    # 2. draw a cube n pixels squared around the sampled pixel
    # 3. calculate the distances and the proportion presence/total
    
    print.noquote("    indexing...")
    # addresses in an array
    address_array <- array(1:(side*side*dim(ch1_t)[3]), 
                           c(side, side, dim(ch1_t)[3]))
    
    # generate the coordinates of pixels in channel1
    ch1_add <- data.frame(which(ch1_t == 1, T))
    colnames(ch1_add) <- c("x", "y", "z")
    # generate the coordinates of pixels in channel2
    ch2_add <- data.frame(which(ch2_t == 1, T))
    colnames(ch2_add) <- c("x", "y", "z")
    
    # avoid pixels too close to the edge
    xbound <- c(size, dim(ch1_t)[1]-(size))
    ybound <- c(size, dim(ch1_t)[2]-(size))
    zbound <- c((size/zstep), dim(ch1_t)[3]-(size/zstep))
    ch1_interior <- ch1_add[ch1_add$x >= xbound[1] & ch1_add$x <= xbound[2] &
                            ch1_add$y >= ybound[1] & ch1_add$y <= ybound[2] &
                            ch1_add$z >= zbound[1] & ch1_add$z <= zbound[2],]
    ch2_interior <- ch2_add[ch2_add$x >= xbound[1] & ch2_add$x <= xbound[2] &
                            ch2_add$y >= ybound[1] & ch2_add$y <= ybound[2] &
                            ch2_add$z >= zbound[1] & ch2_add$z <= zbound[2],]
    
    # determine the sample size of the number of pixels is less than npixels
    if(dim(ch1_interior)[1] < npixels | dim(ch2_interior)[1] < npixels){
      npixels <- min(c(dim(ch1_interior)[1], dim(ch2_interior)[1]))
      print.noquote(paste("no. of pixels <", npixel, ". Reducing npixels to", 
                          npixels, "for this stack"))
    }
    
    # randomly sample the pixels in channel1
    these <- sample(1:dim(ch1_interior)[1], size = npixels)
    # get their addresses
    ch1_pix <- ch1_interior[these,]
    
    # randomly sample the pixels in channel2
    these <- sample(1:dim(ch2_interior)[1], size = npixels)
    # get their addresses
    ch2_pix <- ch2_interior[these,]
    
    # set up a matrix to collect results
    results <- matrix(NA, length(ds), npixels)
    means <- matrix(NA, length(focus)*length(target), length(ds), 
                    dimnames = list(1:(length(focus)*length(target)), ds))
    sds <- means
    
    print.noquote("    calculating...")
    c <- 0
    for(f in focus){
      for(t in target){
        # loop thru to populate the matrix
        for(i in 1:npixels){
          if(f == "ch1"){p <- ch1_pix[i,]}else
            if(f == "ch2"){p <- ch2_pix[i,]}
          # calculate the coordinates of the box
          xrange <- c(p$x-size, p$x+size)
          yrange <- c(p$y-size, p$y+size)
          zrange <- c(p$z-(size/zstep), p$z+(size/zstep))
          # pull the addresses out of the array
          box <- address_array[c(xrange[1]:xrange[2]), 
                               c(yrange[1]:yrange[2]), 
                               c(zrange[1]:zrange[2])]  
          # get the presence absence data for the box
          if(t == "ch1"){id <- ch1_t[box]}else
            if(t == "ch2"){id <- ch2_t[box]}
          # calculate the proportion filled at each distance
          for(l in 1:(length(ds))){
            results[l,i] <- length(which(id[positions[[l]]]==1))/
              length(positions[[l]])
          }
        }
        c <- c+1
        for(i in 1:dim(means)[2]){
          means[c,i] <- mean(results[i,], na.rm = T)
          sds[c,i] <- sd(results[i,], na.rm = T)
        }
        rownames(means)[c] <- paste(f, t, sep = "->")
        rownames(sds)[c] <- paste(f, t, sep = "->")
        print.noquote(paste("        ", f, "->", t))
      }
    }
    
    theseDistances <- cbind(k, melt(means), melt(sds)[,3])
    distance_results <- rbind(distance_results, theseDistances)
  }
  
  # store results in parent directory
  colnames(distance_results) <- c("stack", "direction", "distance", 
                                  "occupancy_mean", "occupancy_sd")
  colnames(productivities) <- c("stack", "strain", "proportion")
  write.csv(distance_results, "distanceOccupancy.csv", row.names = F)
  write.csv(productivities, "productivity.csv", row.names = F)
  
# end function  
}











