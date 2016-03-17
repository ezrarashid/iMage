Tutorial on analysing confocal images microbial from competition experiments
============================================================================

Introduction
------------

This is a brief tutorial on some basic image processing / analysis using
R. The things that I will be showing you here are inspired by Apollo’s
PNAS paper. They include:

1.  Proportion occupancy: For a given pixel in the image, and a
    corresponding distance window, what is the proportion of the
    available space occupied by a given strain? This works by drawing a
    cubed box around the focal pixel, and then within each distance
    window, calculating the proportion of available pixels occupied by a
    given strain. This is something akin to Apollo 2012 (PNAS).

Setting up
----------

The first thing to do is set the working directory in R using `setwd()`.
If you are doing the demo, set your directory to the folder that you
downloaded. I have put the demo data on my computer in
`~/Box Sync/R/iMage001/data/wtAbscess`, but you will need to change the
path to wherever your data is. If you are doing your own images, I am
going to assume for the moment that you have one folder where all of the
.tif files are stored, the .tif files should be separated into
subfolders like in the demo data where each subdirectory contains a set
of files, one for each channel in the image at each slice in the stack.

``` {.r}
path <- "~/Box Sync/R/iMage001/data/wtAbscess" # change this to your path
setwd(path)
```

Note: if you are using windows then replace `/` with `\` and you may
need to type the full path rather than starting with `~`.

You can check that this has worked typing `getwd()`.

``` {.r}
getwd()
```

    ## [1] "/Users/roman/Box Sync/R/iMage001"

Before we get started, you need to load a few libraries into R. If you
are doing this for the first time then you need to install libraries
using the `install.packages()` function e.g. `install.packages("tiff")`.
Once installed, load them like this;

``` {.r}
library(tiff)
library(rtiff)
library(reshape2)
```

    ## Warning: package 'reshape2' was built under R version 3.1.2

``` {.r}
source("~/Box Sync/R/iMage001/iMageFunctions.R")
```

The `source()` command will load the functions that I have defined for
this. Again make sure the path matches the one where you put this
directory. You need to load these two libraries and source the
“iMageFunctions.R” file each time you open R to start this analysis.

Next we need to list all the subdirectories in your `getwd()`, which is
already stored under `path`.

``` {.r}
subdirs <- list.dirs(path, T, F)
subdirs
```

    ## [1] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#1.1"
    ## [2] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#1.2"
    ## [3] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#1.3"
    ## [4] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#2.1"
    ## [5] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#2.2"
    ## [6] "/Users/roman/Box Sync/R/iMage001/data/wtAbscess/wtAbscess#2.3"

This gives R an object over which it can then loop and know which images
are grouped as a stack i.e. grouped into a folder.

Now were ready to go!

1. Proportion occupancy
-----------------------

Were going to first calculate ‘proportion occupancy’ (see intro for the
definition of this). Before we start we need to tell R a few things
about our images.

### 1.1 Setting up some parameters

We need to define a set of parameters that the functions will later use.
For ease I am assuming that all of the images in your `getwd()` are the
same for the parameters defined below. (Basically all stacks the same
size, resolution and z-steps)

``` {.r}
# This stuff tells us about the images
 # (You may need to adjust this for your images)
zstep <- 10        # z-step size in microns
side <- 512        # what is the xy dimension in pixels
pwidth <- 2.4859   # pixel width in microns

# This stuff defines parts of the analysis
 # (for now keep these as they are)
npixel <- 100    # set the number of pixels to be sampled
size <- 30         # set the sampling cube size in pixels
```

Note that the `size` must be a multiple of the `zstep` for the analysis
to run smoothly. You will probably get an error if this is not the case.
Also I am randomly sampling 100 pixels here, see `npixel`. This is a
very small number so that I can test the code quickly. Once you have
everything running, increase this to \~5000 or so to get a more accurate
estimate of the distance occupancy.

### 1.2 Define the channels

Now, in the final part of the setup, we need to tell R which files
belong to which channel. We will do this by entering some text that is
unique to the series of .tif files from that channel. For example, look
at the files in the demo data;

``` {.r}
list.files(subdirs[1])
```

    ##  [1] "wt abscess #1.1 Aa_000.tif" "wt abscess #1.1 Aa_001.tif"
    ##  [3] "wt abscess #1.1 Aa_002.tif" "wt abscess #1.1 Aa_003.tif"
    ##  [5] "wt abscess #1.1 Aa_004.tif" "wt abscess #1.1 Aa_005.tif"
    ##  [7] "wt abscess #1.1 Aa_006.tif" "wt abscess #1.1 Aa_007.tif"
    ##  [9] "wt abscess #1.1 Aa_008.tif" "wt abscess #1.1 Aa_009.tif"
    ## [11] "wt abscess #1.1 Aa_010.tif" "wt abscess #1.1 Aa_011.tif"
    ## [13] "wt abscess #1.1 Aa_012.tif" "wt abscess #1.1 Sg_000.tif"
    ## [15] "wt abscess #1.1 Sg_001.tif" "wt abscess #1.1 Sg_002.tif"
    ## [17] "wt abscess #1.1 Sg_003.tif" "wt abscess #1.1 Sg_004.tif"
    ## [19] "wt abscess #1.1 Sg_005.tif" "wt abscess #1.1 Sg_006.tif"
    ## [21] "wt abscess #1.1 Sg_007.tif" "wt abscess #1.1 Sg_008.tif"
    ## [23] "wt abscess #1.1 Sg_009.tif" "wt abscess #1.1 Sg_010.tif"
    ## [25] "wt abscess #1.1 Sg_011.tif" "wt abscess #1.1 Sg_012.tif"

Looking at these files, we can see that there is a systematic difference
between the two channels in the naming. Files relating to the Aa channel
contain the letters “Aa” and files relating to the Sg channel contain
the letters “Sg”. We’re going to use this to differentiate between the
two channels as follows;

``` {.r}
# unique bit of text pertaining to channel 1
ch1 <- "Aa"

# ... and channel 2
ch2 <- "Sg"
```

### 1.3 Thresholding images

Now we are ready to threshold the images, so that R sees each pixel as a
0 or 1, not an intensity value. For now we are assuming that you have
.tif files with a single colour channel (i.e. monochrome) and that you
have seperated .tif files from your two confocal channels (e.g. strain 1
and 2). In the example data, the files are monochrome and they are
already binary (only contain values of 0 and 1). In this case, we don’t
need to threshold them but we still need a processing step that will
convert the .tiff files into an easily manipulated R file (containing an
Array, or the whole z-stack). I have written a function for this called
`tiffToArray()`, which we will use as follows.

``` {.r}
tiffToArray(subdirs = subdirs, 
            side = side, 
            channels = c(ch1, ch2), 
            threshold = F
            )
```

The first argument `subdir` is the list of subdirectories where we want
it to look (we designated this earlier). The second argument `side` is
the number of pixels on each side of each .tiff. Again we designated
this in 1.1. The third argument is called `channels` and takes a vector
of values denoting the unique bit of text in the filenames for each
channel that we identified in 1.2.

When you run this command, you should see a whole bunch of output in
your terminal and some new files in your parent directory. As I
mentioned, the demo dataset is already thresholded but in case you need
to threshold your images, just replace the argument `threshold = F` with
`threshold = T`.

### 1.4 Calculating proportion occupancy

The last thing to do before running the calculations is to tell R where
the new thresholded R files are. Use the function `list.files()` as
follows;

``` {.r}
thresholdedStackFiles <- list.files(path, "thresholdedStack", full.names = T)
```

Now we’re ready to run the number cruncher function which will run
through all your stacks and calculate proportion occupancy (inspired by,
but not identical to Stacey 2012). Briefly the process is as follows:

1.  A sampling box is constructed (according to parameter `size` that
    you set)
2.  A z-stack (actually 2, one for each channel) is loaded into memory.
3.  For `npixel` randomly selected pixels the box is centred on that
    pixel.
4.  For each placement of the sampling box, the algorithm moves from the
    centre to the outer corners in 2uM steps.
5.  Within each 2uM window, we record the proportion of the total space
    available that is occupied by the target strain. (e.g. a value of 1
    would mean that all of the available space is occupied by the strain
    of interest.

The most important bit is what happens in (e). This is similar, but not
identical to what’s measured in Stacey 2012.

Run the function on the demo dataset as follows:

``` {.r}
proportionOccupancy(files = thresholdedStackFiles, 
                    channels = c(ch1, ch2),
                    zstep = zstep,
                    side = side,
                    pwidth = pwidth,
                    npixel = npixel,
                    size = size,
                    focus = c("ch1", "ch2"),
                    target = c("ch1", "ch2"),
                    cleanup = T
                    )
```

The argument `files` is telling R where the newly created files are
(created using `tiffToArray()`). The argument `channels` is passing
information we defined previously about some unique filename feature.
The arguments `zstep`, `side`, `pwidth`, `npixel`, and `size` all
contain information that we defined earlier in 1.1 and some of which was
also used in `tiffToArray()`.

The two arguments `focus`, and `target` determine which channels we
would like to consider as the ‘focal’ and ‘target’ strains. For example
we might be interested only in the distance of channel 2 from a focus of
pixels in channel 1. The default is to do the calculations for all 4
combinations. You don’t need to include these arguments if that’s what
you want.

There is some built in default behaviour for example if you have
specified a value of `npixel` and there are fewer than `npixel` in a
particular channel then the program will look at all of the pixels once,
so your sample size for that stack will be a little lower. You will get
a warning in the output if this happens. Generally speaking it should
not be a problem if you are replicating at the level of the z-stack.

The argument `cleanup` determines whether the files created by
`tiffToArray` are deleted.

When you run this on the demo data, you should see another load of
output as it runs. When it stops, you should have two new files
“distanceOccupancy.csv” and “productivity.csv” in your directory. Thats
our data, lets plot it.

### 1.5 Plotting the results

First read in the new data files.

``` {.r}
prod <- read.csv("productivity.csv")
dist <- read.csv("distanceOccupancy.csv")
```

Lets first plot the productivities. This is just the proportion of total
available space in the stack that is occupied by a particular strain. A
value of 1 would mean that strain had occupied the entire stack (every
pixel). First we use the function `subset()` to create two subsets of
the data depending on which strain it is. Then we use the `plot()`
function to plot the points. What results is a scatter plot of the
productivities. Each data point represents one confocal z-stack.

``` {.r}
aa_prod <- subset(prod, strain == "Aa")
sg_prod <- subset(prod, strain == "Sg")
plot(aa_prod$proportion ~ sg_prod$proportion,
     ylab = "Aa productivity",
     xlab = "Sg productivity")
```

Now lets plot the distance occupancy profiles. If you look at the `dist`
object that we have just read into memory. Clue: useful functions for
looking at R objects are `head()` and `str()`. For example try
`head(dist)` now. You’ll see that there is a variable called direction.
This tells us in which direction we are looking e.g. from a ‘focal’
pixel to the surrounding ‘target’ pixels. At the moment this is talking
about channels but we want to know about strains. I’m using the `gsub()`
function to substitute the text about channels for some other text about
strains as follows:

``` {.r}
dist$direction <- gsub("ch1", "Aa", dist$direction)
dist$direction <- gsub("ch2", "Sg", dist$direction)
```

Now lets make a plot. I am using the `par()` function to split the
plotting space into 4 sections. Then I am defining a variable
`directions` which will store the different combinations of focal and
target strain in our data. Then I am using a loop to loop thru
`directions` and make a plot for each one. In each plot I am subsetting
the data again and then making an empty plot first with
`plot(type = "n")`. The for each stack I am subsetting the data yet
again and then adding a line with `lines()` for each z-stack.

``` {.r}
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
```

Admittedly this is quite a long winded way of plotting, if you have
already started learning something better like ggplot2 then stick with
that. The plots look quite ugly but we can easily make them more
beautiful another time.

What do you notice? There is quite a lot of enrichment of Aa close to Aa
and Sg close to Sg (first and last panels) but not the other way. This
indicates that strains are clustered, not so surprising given fission.
Anyway we can leave the biology for another day but hopefully some
analysis along these lines will help us go quickly from images through
to interpretable output in a reasonably automated (if a little slow)
way.

:-)

Just to recap, your final script could look something like this:

``` {.r}
path <- "~/Box Sync/R/iMage001/data/wtAbscess"
setwd(path)

library(tiff)
library(rtiff)
library(reshape2)
source("~/Box Sync/R/iMage001/iMageFunctions.R")

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
thresholdedStackFiles <- list.files(path, "thresholdedStack", full.names = T)

# do the number crunching on proportion occupancy
proportionOccupancy(files = thresholdedStackFiles, 
                    channels = c(ch1, ch2),
                    zstep = zstep,
                    side = side,
                    pwidth = pwidth,
                    npixel = npixel,
                    size = size,
                    focus = c("ch1", "ch2"),
                    target = c("ch1", "ch2"),
                    cleanup = T
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
```
