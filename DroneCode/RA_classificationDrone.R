# Red alder classification project with UW using drone lidar data
#

# Introductory comments ---------------------------------------------------
# This is the updated code from the original RA classification project. The code
# has been updated to use drone lidar data instead of airborne data. The drone
# data are very high density (1000+ returns/square meter). However, the intensity
# values are noisy and may be less useful than those from the airborne data.
#
# There are 2 files for the trees for each plot. The first has the "turning point" tree (this is
# the tree used as the origin for tree locations) and the second has the bulk of the
# trees on the plot. Plots are square. The data tables are different for the two
# files with the TP file having fewer columns. I only checked 2 files but the format
# looks consistent between the 2 files. I don't really know what the fields represent
# but I don't want to lose anything so I will add dummy fields to the TP trees. There
# are also some problems caused by limitations with shapefile column names. Not sure why
# there appear to be duplicate columns but there are a few rows that have different
# values in the "duplicate" columns.
#
# Lidar data for this work is organized by plot. FOr the previous project, all data were
# in a single folder. For this project, the calls to functions are a little more difficult
# since we have to pass the folder name for the data.
#
# equation for sample circle from Ally:
# circle radius = r = 0.28t - 6.62
# Where r = radius and t = tree height (all in feet). For this we actually don't have tree
# heights so we will use the highest return height.
#
# Units for drone lidar are meters for XY and Z so the equation above needs some work. It is
# easiest to convert height and radius going into and out of the equation.
#
# Code needs other modifications to work with meters.
#
# for processPlot, changed assumedHeight = 40 and minRadius = 2.0 to assumedHeight = 12.19 and minRadius = 0.61


# mergeTrees function -----------------------------------------------------
mergeTrees <- function (
  plotName = NULL,
  folder = NULL,
  test = TRUE,
  logCon = NA
) {
  # check parameters
  if (is.null(plotName) || is.null(folder)) {
    stop("You must specify plotName and folder")
  }
  
  # with test = TRUE, the code will catch field errors related to missing Anomaly_Nu field and "fix" them
  # with test = FALSE, the code will fail when the Anomaly_Nu field is missing
  
  # read the tree files...trees and TP
  trees <- st_read(paste0(folder, plotName, "_Trees_Project.shp"))
  TPtrees <- st_read(paste0(folder, plotName, "_TP_Project.shp"))
  
  if (test) {
    newtrees <- tryCatch(trees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text", "Anomaly_Nu")], error = function(e) {NA})
    if (!is.object(newtrees)) {
      newtrees <- trees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text")]
      newtrees$Anomaly_Nu <- 0
      
      if (is.object(logCon)) {
        writeLines("  Missing Anomaly_Nu field in Trees", logCon)
      }
    }
  
    newTPtrees <- tryCatch(TPtrees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text", "Anomaly_Nu")], error = function(e) {NA})
    if (!is.object(newTPtrees)) {
      newTPtrees <- TPtrees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text")]
      newTPtrees$Anomaly_Nu <- 0
    
        if (is.object(logCon)) {
        writeLines("  Missing Anomaly_Nu field in TP", logCon)
      }
    }
  } else {  
    newtrees <- trees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text", "Anomaly_Nu")]
    newTPtrees <- TPtrees[, c("Tag", "Species", "OldTag", "DBH_cm", "Note_num", "Note_text", "Anomaly_Nu")]
  }
  
  t <- rbind(newtrees, newTPtrees)
  
  # reproject
  t <- st_transform(t, crs = st_crs("EPSG:26910"))
  
  rm(trees)
  rm(TPtrees)
  rm(newtrees)
  rm(newTPtrees)
  
  return(t)
}

# clipCrownCircle function ------------------------------------------------
clipCrownCircle <- function (
  plot,
  treeID,
  X,
  Y,
  radius,
  lastReturns = FALSE,
  inputFiles = NULL,
  outputFolder = NULL,
  groundFile = NULL
) {
  # check parameters
  if (is.null(inputFiles) || is.null(outputFolder) || is.null(groundFile)) {
    stop("You must specify inputFiles, outputFolder, and groundFile")
  }
  
  # create the output folder
  dir.create(outputFolder, showWarnings = FALSE)

  lastFlag <- ""
  if (lastReturns) lastFlag <- "/return:L /zero"
  
  cmd <- paste("clipdata64"
              , paste0("\"/ground:", groundFile, "\"")
              , paste0("/shape:1")
              , "/height"
              , lastFlag
              , paste0("\"", inputFiles, "\"")
              , paste0("\"", outputFolder, plot, "_Tree_", treeID, ".las\"")
              , X - radius
              , Y - radius
              , X + radius
              , Y + radius
  )
  
#  cat(cmd, "\n")
  shell(cmd)
}

# computeMetrics function -------------------------------------------------
computeMetrics <- function (
  plot,
  heightThreshold = 1.37,
  coverThreshold = 1.37,
  firstReturns = TRUE,
  inputFiles = NULL,
  outputFolder = NULL,
  outputFileBaseName = "temp_metrics.csv",
  deleteClips = FALSE
) {
  # check parameters
  if (is.null(inputFiles) || is.null(outputFolder)) {
    stop("You must specify inputFiles and outputFolder")
  }
  
  first <- ""
  if (firstReturns) first <- "/first"
  
  cmd <- paste("cloudmetrics"
               , "/new"
               , "/rid"
               , paste0("/above:", coverThreshold)
               , paste0("/minht:", heightThreshold)
               , first
               , paste0("\"", inputFiles, "\"")
               , paste0("\"", outputFolder, plot, "_", outputFileBaseName, "\"")
  )
  
  if (shell(cmd) == 0) {
    if (deleteClips) {
      unlink(inputFiles)
    }
  }
}

# processPlot function ----------------------------------------------------
processPlot <- function(
  plotID,
  assumedHeight = 12.19,
  minRadius = 0.61,
  radiusFudgeFactor = 1.0,
  logCon = NULL,
  treeDataFolder = NULL,
  inputFiles = NULL,
  outputFolder = NULL,
  groundFile = NULL,
  tempMetricsFileBase = "temp_metrics.csv"
) {
  # check parameters
  if (is.null(inputFiles) || is.null(outputFolder) || is.null(treeDataFolder) || is.null(groundFile)) {
    stop("You must specify inputFiles, outputFolder, treeDataFolder, and groundFile")
  }
  
  #  plotID <- "B1E3_TP_25_100"     # testing only
  
  if (!is.null(logCon)) {
    writeLines(paste("Processing plot:", plotID), logCon)
  }
  
  trees <- mergeTrees(plotID, treeDataFolder, logCon = logCon)
  
#  mapview(trees, zcol = "Species")
  
  # assume a height of 40 feet to create an initial sample...must be larger than ~24 feet 
  # to produce a radius > 0
  trees$FalseHeight <- assumedHeight
  #trees$FalseRadius <- 0.28 * trees$FalseHeight - 6.62
  trees$FalseRadius <- (0.28 * (trees$FalseHeight * 3.2808) - 6.62) / 3.2808
  
  # add coordinates as separate columns...there is a better way to do this but I 
  # don't know it...
  coords <- st_coordinates(trees)
  trees <- cbind(trees, coords)
  
  # clip the sample for each tree...res is a list of error codes
  # expect them all to be 0 it clips worked
  res <- apply(trees, MARGIN = 1, FUN = function(x) {
    clipCrownCircle(plotID
                    , x$Tag
                    , x$X
                    , x$Y
                    , x$FalseRadius
                    , inputFiles = inputFiles
                    , outputFolder = outputFolder
                    , groundFile = groundFile
    )
  })
  
  # compute metrics
  res <- computeMetrics(plotID, deleteClips = TRUE, inputFiles = paste0(outputFolder, "*.las"), outputFolder = outputFolder)
  
  # read metrics and extract highest point in the sample
  metrics <- read.csv(paste0(outputFolder, plotID, "_", tempMetricsFileBase), stringsAsFactors = FALSE)
  
  # compute new circle radius using P99 height
  trees <- merge(trees, metrics[, c(1, 50)], by.x = "Tag", by.y = "Identifier")
  trees$LidarHt <- trees$Elev.P99
  trees$TrueRadius <- (0.28 * (trees$LidarHt * 3.2808) - 6.62) / 3.2808 * radiusFudgeFactor
  
  # condition the radius...radius equation goes negative for heights < 23.64 feet (7.2 m)
  trees$TrueRadius[trees$TrueRadius < minRadius] <- minRadius
  
  # delete the Elev.P99 column...the value is in the LidarHt column. If left in the data frame,
  # we end up with 2 P99 values with ".x" and ".y" modifiers
  trees$Elev.P99 <- NULL
  
  # clip new samples
  res <- apply(trees, MARGIN = 1, FUN = function(x) {
    clipCrownCircle(plotID
                    , x$Tag
                    , x$X
                    , x$Y
                    , x$TrueRadius
                    , inputFiles = inputFiles
                    , outputFolder = outputFolder
                    , groundFile = groundFile
    )
  })
  
  # compute metrics...first returns only
  res <- computeMetrics(plotID, deleteClips = FALSE, inputFiles = paste0(outputFolder, "*.las"), outputFolder = outputFolder)
  
  # append to features
  # read metrics, change column labels and append to tree records
  metrics <- read.csv(paste0(outputFolder, plotID, "_", tempMetricsFileBase), stringsAsFactors = FALSE)
  colnames(metrics) <- paste0("First.", colnames(metrics))
  treeMetrics <- merge(trees, metrics[, -c(2:3)], by.x = "Tag", by.y = "First.Identifier")
  
  # compute metrics...all returns
  res <- computeMetrics(plotID, firstReturns = FALSE, deleteClips = TRUE, inputFiles = paste0(outputFolder, "*.las"), outputFolder = outputFolder)
  
  # append to features
  # read metrics and append to tree records
  metrics <- read.csv(paste0(outputFolder, plotID, "_", tempMetricsFileBase), stringsAsFactors = FALSE)
  treeMetrics <- merge(treeMetrics, metrics[, -c(2:3)], by.x = "Tag", by.y = "Identifier")

  # clip last returns
  res <- apply(trees, MARGIN = 1, FUN = function(x) {
    clipCrownCircle(plotID
                    , x$Tag
                    , x$X
                    , x$Y
                    , x$TrueRadius
                    , lastReturns = TRUE
                    , inputFiles = inputFiles
                    , outputFolder = outputFolder
                    , groundFile = groundFile
    )
  })
  
  # compute metrics using last returns
  res <- computeMetrics(plotID, heightThreshold = 0.0, firstReturns = FALSE, deleteClips = TRUE, inputFiles = paste0(outputFolder, "*.las"), outputFolder = outputFolder)
  
  # append to features
  # read metrics and append to tree records
  metrics <- read.csv(paste0(outputFolder, plotID, "_", tempMetricsFileBase), stringsAsFactors = FALSE)
  
  # only want a few metrics...P05, P10, P50, mean
  metrics <- metrics[, c("Identifier", "Elev.P01", "Elev.P05", "Elev.P10", "Elev.P50", "Elev.mean")]
  
  # change names
  colnames(metrics) <- c("Identifier", "Last.Elev.P01", "Last.Elev.P05", "Last.Elev.P10", "Last.Elev.P50", "Last.Elev.mean")
  treeMetrics <- merge(treeMetrics, metrics, by.x = "Tag", by.y = "Identifier")

  # replace "." in column names with "_" and write trees as geopackage
  # Arc doesn't like the "." but it could be geopackages in general
  names(trees) <- gsub("\\.", "_", names(trees))
  st_write(trees, paste0(plotID, "_Trees.gpkg"), "Trees", delete_dsn = TRUE, delete_layer = TRUE)
  
  # buffer trees using the final circle radius and save polygons
  treePolys <- st_buffer(trees, dist = trees$TrueRadius)
  
  # moved this to write tree circles after adding the overlap info
#  st_write(treePolys, paste0(plotID, "_Trees.gpkg"), "SampleCircles", delete_dsn = FALSE, delete_layer = TRUE)
  
#  mapview(treePolys, zcol = "Species")
  
  # delete temp metrics file
  unlink(paste0(outputFolder, plotID, "_", tempMetricsFileBase))
  
  # figure out which tree circles overlap other circles
  t <- st_overlaps(treePolys, sparse = FALSE)
  treePolys$Overlapped <- (colSums(t) > 0)
  
#  mapview(treePolys, zcol = "Overlapped")
  
  # figure out if overlapping circles are the same species...sparse = TRUE gives lists of overlapping circles (if any)
  t <- st_overlaps(treePolys, sparse = TRUE)
  
  # loop seems like the easiest way to do this but I'm sure there is a "vectorized" way to do it
  treePolys$OverlapSpeciesDifferent <- FALSE
  for (i in 1:nrow(treePolys)) {
    # compare species for overlapping circles and flag if any are different
    if (length(t[[i]]) > 0) {
      for (j in 1:length(t[[i]])) {
        if (treePolys$Species[i] != treePolys$Species[t[[i]][j]]) {
          treePolys$OverlapSpeciesDifferent[i] <- TRUE
          break
        }
      }  
    }
  }
  #treePolys$OverlapSameSpecies <- (colSums(t) > 0)
  
#  mapview(treePolys, zcol = "OverlapSpeciesDifferent")

  st_write(treePolys, paste0(plotID, "_Trees.gpkg"), "SampleCircles", delete_dsn = FALSE, delete_layer = TRUE)
  
  # add overlap info to metrics
  treeMetrics$Overlapped <- treePolys$Overlapped
  treeMetrics$OverlapSpeciesDifferent <- treePolys$OverlapSpeciesDifferent
  
  # drop geometry...causes problems when reading into excel
  treeMetrics <- st_set_geometry(treeMetrics, NULL)
  
  # rearrange columns
  treeMetrics <- treeMetrics[, c(1:13, 221:222, 14:220)]
  
  # write metrics as CSV
  write.csv(treeMetrics, paste0(plotID, "_TreeMetrics.csv"), row.names = FALSE)
}




# Main code section -------------------------------------------------------
# *****************************************************************************
# Basic setup stuff
# *****************************************************************************
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(randomForest)
library(caret)
library(yaImpute)
library(mapview)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)


setwd("E:/Backup/R_Stuff/AlderClassification/DroneResults")

# I grabbed all of the plot data from the shared google drive 8/11/2021
# there are 2 shapefiles for each plot. One has the turning point and the other has the trees
# the plots.txt file lists the identifiers for the plots
dataFolder <- "E:/Backup/R_Stuff/RAClassification/Georeferenced_tree_plots_2020_sampling/"

# this is the direct link to the shared google drive
#dataFolder <- "I:/.shortcut-targets-by-id/1j0spS2MzKbPSsif_DEv28S0AjAenk_Hj/Georeferenced_tree_plots_2020_sampling/"

# this is a folder with a few plots for testing
#dataFolder <- "g:/R_Stuff/RAClassification/data/"

# folder for temporary work...could be anywhere...not sure why I put it where it is
#tempClipFolder <- "E:/Backup/temp/Keven_polys/tempclips/"
#tempMetricsFileBase <- "temp_metrics.csv"

# this folder has metrics for all and first returns...different layer names
# on my RESULTS013 drive
rasterFolder <- "J:/ONRC_SOLDUC/AP_run_5ft/Products_Sappho_2021-04-14/FINAL_Sappho_2021-04-14/Metrics_5FEET/"

# this folder has metrics for last returns only
lastRasterFolder <- "J:/ONRC_SOLDUC/AP_run_5ft_LAST/Products_Sappho_2021-04-26/FINAL_Sappho_2021-04-26/Metrics_5FEET/"

resultsFolder <- "g:/R_Stuff/RAClassification/DroneResults/"

#minRadius <- 2.0 # feet
minRadius <- 0.61 # meters
radiusFudgeFactor <- 1.0 # multiplier...no units

# *****************************************************************************
# process plots
# *****************************************************************************
# test plot
#processPlot("B1E3_TP_25_100")
#TPtrees <- st_read(paste0(dataFolder, "B4E3_TP_50_0", "_TP_Project.shp"))

con <- file("../plots.txt")
#con <- file("shortlist_plots.txt")
plotList <- readLines(con)
close(con)

# read lidar plot IDs from excel workbook...only has one sheet
t <- read_excel("../PlotID.xlsx")

# build plot identifier so we can match to lidar data folders
plotIDs <- substr(plotList, 1, 4)

lidarIDs <- t[t$LTEPPlotID %in% plotIDs, ]

# clear the status file
if (file.exists("PlotProcessing.txt")) {
  unlink("PlotProcessing.txt")
}

logCon <- file("PlotProcessing.txt", "wt")

for (i in 1:length(plotList)) {
  # get the 4 letter plot identifier
  plot <- substr(plotList[i], 1, 4)
  
  # see if we have lidar data
  lidarID <- lidarIDs$DronePlotID[which(lidarIDs$LTEPPlotID == plot)]
  
  processPlot(plotList[i], assumedHeight = 37 / 3.2808, logCon = logCon,
              treeDataFolder = dataFolder,
              inputFiles = paste0("E:/2022_DroneLidar/Sappho LTEP/Plot ", lidarID, "/*.laz"),
              outputFolder = "C:/Temp/treeclips/",
              groundFile = paste0("E:/2022_DroneLidar/Sappho LTEP/Plot ", lidarID, "/ground/ground.dtm")
  )
}

close(logCon)

#processPlot("B3E2_TP_75_50", assumedHeight = 37)

# *****************************************************************************
# merge all circles for individual plots
# *****************************************************************************
con <- file("plots.txt")
plotList <- readLines(con)
close(con)

for (i in 1:length(plotList)) {
  if (file.exists(paste0(plotList[i], "_Trees.gpkg"))) {
    # read plot file
    tm <- st_read(paste0(plotList[i], "_Trees.gpkg"), "SampleCircles")
    
    # add a plot identifier and reorder
    tm$PlotID <- plotList[i]
    tm <- tm[, c(ncol(tm), 1:(ncol(tm) - 1))]
    
    if (i == 1) {
      allTrees <- tm
    } else {
      allTrees <- rbind(allTrees, tm)
    }
  }
}
st_write(allTrees, "AllPlots_Trees.gpkg", "SampleCircles", delete_dsn = TRUE, delete_layer = TRUE)
all_Trees <- st_read("AllPlots_Trees.gpkg")

# filter based on overlap information...do this later
#modelTrees <- dplyr::filter(allTrees, !Overlapped | !OverlapSpeciesDifferent)
modelTrees <- allTrees

#mapview(allTrees, zcol = "Species")
#mapview(modelTrees, zcol = "Species")

# *****************************************************************************
# merge all TreeMetrics files for individual plots and extract non-overlapped 
# trees and trees overlapped by same species trees, compute new metrics
# and write CSV file with all trees and metrics
# *****************************************************************************
con <- file("plots.txt")
plotList <- readLines(con)
close(con)

for (i in 1:length(plotList)) {
  if (file.exists(paste0(plotList[i], "_TreeMetrics.csv"))) {
    # read plot file
    tm <- read.csv(paste0(plotList[i], "_TreeMetrics.csv"), stringsAsFactors = FALSE)
    
    # add a plot identifier and reorder
    tm$PlotID <- plotList[i]
    tm <- tm[, c(ncol(tm), 1:(ncol(tm) - 1))]
    
    if (i == 1) {
      allMetrics <- tm
    } else {
      allMetrics <- rbind(allMetrics, tm)
    }
  }
}

# filter based on overlap information...do this later
#modelMetrics <- dplyr::filter(allMetrics, !Overlapped | !OverlapSpeciesDifferent)
modelMetrics <- allMetrics

# add relative height percentiles...relative to P95...for all return metrics
modelMetrics$Elev.RP01 <- modelMetrics$Elev.P01 / modelMetrics$Elev.P95
modelMetrics$Elev.RP05 <- modelMetrics$Elev.P05 / modelMetrics$Elev.P95
modelMetrics$Elev.RP10 <- modelMetrics$Elev.P10 / modelMetrics$Elev.P95
modelMetrics$Elev.RP20 <- modelMetrics$Elev.P20 / modelMetrics$Elev.P95
modelMetrics$Elev.RP25 <- modelMetrics$Elev.P25 / modelMetrics$Elev.P95
modelMetrics$Elev.RP30 <- modelMetrics$Elev.P30 / modelMetrics$Elev.P95
modelMetrics$Elev.RP40 <- modelMetrics$Elev.P40 / modelMetrics$Elev.P95
modelMetrics$Elev.RP50 <- modelMetrics$Elev.P50 / modelMetrics$Elev.P95
modelMetrics$Elev.RP60 <- modelMetrics$Elev.P60 / modelMetrics$Elev.P95
modelMetrics$Elev.RP70 <- modelMetrics$Elev.P70 / modelMetrics$Elev.P95
modelMetrics$Elev.RP75 <- modelMetrics$Elev.P75 / modelMetrics$Elev.P95
modelMetrics$Elev.RP80 <- modelMetrics$Elev.P80 / modelMetrics$Elev.P95
modelMetrics$Elev.RP90 <- modelMetrics$Elev.P90 / modelMetrics$Elev.P95

# add relative intensity percentiles...relative to P95...for first returns
modelMetrics$First.Int.RP01 <- modelMetrics$First.Int.P01 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP05 <- modelMetrics$First.Int.P05 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP10 <- modelMetrics$First.Int.P10 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP20 <- modelMetrics$First.Int.P20 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP25 <- modelMetrics$First.Int.P25 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP30 <- modelMetrics$First.Int.P30 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP40 <- modelMetrics$First.Int.P40 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP50 <- modelMetrics$First.Int.P50 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP60 <- modelMetrics$First.Int.P60 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP70 <- modelMetrics$First.Int.P70 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP75 <- modelMetrics$First.Int.P75 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP80 <- modelMetrics$First.Int.P80 / modelMetrics$First.Int.P95
modelMetrics$First.Int.RP90 <- modelMetrics$First.Int.P90 / modelMetrics$First.Int.P95

# add last/first return metrics
modelMetrics$Last.P05.over.First.P95 <- modelMetrics$Last.Elev.P05 / modelMetrics$First.Elev.P95
modelMetrics$Last.P10.over.First.P90 <- modelMetrics$Last.Elev.P10 / modelMetrics$First.Elev.P90
modelMetrics$Last.mean.over.First.mean <- modelMetrics$Last.Elev.mean / modelMetrics$First.Elev.mean
modelMetrics$First.P95.minus.Last.P05 <- modelMetrics$First.Elev.P95 - modelMetrics$Last.Elev.P05
modelMetrics$First.P90.minus.Last.P10 <- modelMetrics$First.Elev.P90 - modelMetrics$Last.Elev.P10
modelMetrics$First.mean.minus.Last.mean <- modelMetrics$First.Elev.mean - modelMetrics$Last.Elev.mean

# add new predictor variables...all returns
modelMetrics$Elev.RP25_75 <- modelMetrics$Elev.RP25 / modelMetrics$Elev.RP75
modelMetrics$Elev.RP10_90 <- modelMetrics$Elev.RP10 / modelMetrics$Elev.RP90

# write data for model development
write.csv(modelMetrics, paste0("AllPlots", "_ModelMetrics.csv"), row.names = FALSE)

# summary info...circle radius by species
summary(modelMetrics$TrueRadius)
summary(modelMetrics$TrueRadius[modelMetrics$Species == "ALRU"])
summary(modelMetrics$TrueRadius[modelMetrics$Species == "PSME"])
summary(modelMetrics$TrueRadius[modelMetrics$Species == "TSHE"])
summary(modelMetrics$TrueRadius[modelMetrics$Species == "PISI"])
summary(modelMetrics$TrueRadius[modelMetrics$Species == "ACCI"])
summary(modelMetrics$TrueRadius[modelMetrics$Species == "RHPU"])

table(modelMetrics$Species)


# Modeling with RF --------------------------------------------------------


# *****************************************************************************
# fit a classification model...OK to jump directly to this code after setup block
# once the CSV file is generated
# *****************************************************************************
# this variable controls some of the code below that sets up modeling
modelType <- "ALRU"
#modeType <- "PSME"

seed <- 34567
set.seed(seed)

# read data
trainData <- read.csv(paste0("AllPlots", "_ModelMetrics.csv"), stringsAsFactors = FALSE)

# make species a factor...necessary for random forest to do classification
trainData$Species <- as.factor(trainData$Species)

# drop trees with anomalies:
# can drop all trees with Anomaly_Nu of 0
# OR drop trees with values of 1 or 2
# Anomaly_Nu:
# 0: No problem
# 1: Dead tree
# 2: Leaning tree
# 3: Tree shares base with other tree/forked tree
# 4: Tree is on stump
# 5: Tree is less than 0.5 meters away from another tree
# 6: Check in field (location doesn't make sense for some reason or another)
#trainData <- dplyr::filter(trainData, Anomaly_Nu == 0)
trainData <- dplyr::filter(trainData, Anomaly_Nu == 0 | Anomaly_Nu == 3 | Anomaly_Nu == 4 | Anomaly_Nu == 5)

# drop trees with DBH < 10cm
trainData <- dplyr::filter(trainData, DBH_cm >= 10)

# filter using overlap information
# *********** used for initial testing but found that including the overlapped trees
# gave better predictions for the hold-out testing data
#trainData <- dplyr::filter(trainData, !Overlapped | !OverlapSpeciesDifferent)

# get counts by Species
table(trainData$Species)

# show some summary info for variables used for RF model to classify ALRU
trainData %>%
  group_by(Species) %>%
  summarize(
    Count = n(),
    AveP90 = mean(Elev.P90, na.rm = TRUE),
    FirstIntP80 = mean(First.Int.P80, na.rm = TRUE),
    IntP20 = mean(Int.P20, na.rm = TRUE),
    IntIQ = mean(Int.IQ, na.rm = TRUE),
    FirstLastMeanDiff = mean(First.mean.minus.Last.mean, na.rm = TRUE)
  )

# species after dropping trees with anomalies and < 10cm DBH
#ACCI ALRU PISI PSME RHPU TSHE  UNK UNKN 
#0    150    8  475    1  356    0    0 

# species counts after dropping trees overlapped by different species
#ACCI ALRU PISI PSME RHPU TSHE  UNK UNKN 
#0   73    4  306    1  268    0    0 

if (modelType == "ALRU") {
  # create a new "species" code that groups ACCI, PISI, and RHPU into OTHER
  # also lumped PSME and TSHE into CONIFER type
  trainData$SpeciesGroup[trainData$Species == "ALRU"] <- "ALRU"
  trainData$SpeciesGroup[trainData$Species == "PSME"] <- "CONIFER"
  trainData$SpeciesGroup[trainData$Species == "TSHE"] <- "CONIFER"
  trainData$SpeciesGroup[trainData$Species == "PISI"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "ACCI"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "RHPU"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "UNK"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "UNKN"] <- "OTHER"
} else {
  # reclass for PSME/TSHE...I tested this and the model has accuracy = ~79%
  # all of the important metrics are height related
  trainData$SpeciesGroup[trainData$Species == "ALRU"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "PSME"] <- "PSME"
  trainData$SpeciesGroup[trainData$Species == "TSHE"] <- "TSHE"
  trainData$SpeciesGroup[trainData$Species == "PISI"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "ACCI"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "RHPU"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "UNK"] <- "OTHER"
  trainData$SpeciesGroup[trainData$Species == "UNKN"] <- "OTHER"
}

# show some summary info
trainData %>%
  group_by(Species) %>%
  summarize(
    count = n(),
    aveP70 = mean(Int.P70, na.rm = TRUE),
    aveCovAboveMean = mean(Percentage.all.returns.above.mean, na.rm = TRUE)
    )

# drop OTHER group...only has a few individuals with DBH > 10cm
trainData <- dplyr::filter(trainData, SpeciesGroup != "OTHER")

trainData$SpeciesGroup <- as.factor(trainData$SpeciesGroup)

# get counts by SpeciesGroup...with DBH >= 10cm, we only have 9 OTHER
table(trainData$Species)
table(trainData$SpeciesGroup)

summary(trainData$TrueRadius)

# drop some metrics: return counts
# column 3 has the actual species, column 258 has the species group identifier
trainData <- trainData[, -c(1:30, 103:104, 118:131, 204:205, 212:217)]   # by SpeciesGroup

# drop records for trees with no returns...RP## = NA
# this drops 3 trees...1 had 1 point and the other 2 0 points
trainData <- na.omit(trainData)

# show some summary info
trainData %>%
  group_by(SpeciesGroup) %>%
  summarize(
    count = n(),
    aveP90 = mean(Elev.P90, na.rm = TRUE),
    aveCovAboveMean = mean(Percentage.all.returns.above.mean, na.rm = TRUE)
  )

#plot(trainData$Elev.RP25 / trainData$Elev.RP75, trainData$Percentage.all.returns.above.mean)
#plot(trainData$Elev.RP25 / trainData$Elev.RP75, trainData$SpeciesGroup)
#plot(trainData$Elev.RP10 / trainData$Elev.RP90, trainData$Percentage.all.returns.above.mean)
#plot(trainData$Elev.RP10 / trainData$Elev.RP90, trainData$SpeciesGroup)

# *****************************************************************************
# *****************************************************************************
# *****************************************************************************
# *****************************************************************************
# *****************************************************************************
# at this point, we have good data


# *****************************************************************************
# most of the analysis done for the paper was moved to the ComparisonTable.R
# file. The mapping code is still below. One flaw in my logic is that I ran the
# RF model 30 times with different data splits and random number seeds. This
# gives very stable estimates of model performance but I can't use the 30 runs
# to "apply" the model to map alder
# *****************************************************************************

# for PSME/TSHE classification, may want/need to drop most of the metrics related to elevation...
# PSME is taller, on average, so we expect to be able to classify based on height metrics alone
#SpeciesGroup     count   aveP90 aveCovAboveMean
#<fct>            <int>   <dbl>           <dbl>
#  1 PSME           474   36.4            46.7
#  2 TSHE           354   25.5            41.8
plot(trainData$SpeciesGroup, trainData$Elev.P99, xlab = "Species group", ylab = "P99 height (ft)")

#plot(trainData$SpeciesGroup, trainData$Int.IQ, xlab = "Species group")

par(mfrow = c(2, 2), oma = c(1, 1, 0, 0), mar = c(2, 4.5, 1, 1))
plot(trainData$SpeciesGroup, trainData$First.Int.P80, xlab = "Species group", ylab = expression(80^th~percentile~of~first~return~intensity))
plot(trainData$SpeciesGroup, trainData$Int.P20, xlab = "Species group", ylab = expression(20^th~percentile~of~all~return~intensity))
plot(trainData$SpeciesGroup, trainData$Int.IQ, xlab = "Species group", ylab = "All return intensity IQ distance")
plot(trainData$SpeciesGroup, trainData$First.mean.minus.Last.mean, xlab = "Species group", ylab = "First return mean height - last return mean height (ft)")
par(mfrow = c(1, 1))

# work on improving plots...still having problems with alignment due to superscripts in Y-axis labels
# removed the superscripted parts
fs <- 14
ff <- "sans"
P80 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=First.Int.P80, fill=SpeciesGroup)) + 
  geom_boxplot() +
  geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "First.Int.P80"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "First.Int.P80"])) / 2, linetype="dashed", color = "black") +
  #  geom_violin(draw_quantiles = c(.50), trim = FALSE) +
  scale_x_discrete(labels = NULL) +
  ylab("80th percentile of first return intensity") +
  #  ylab(expression(80th~percentile~of~first~return~intensity)) +
  xlab("") +
  theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.5,0,0.1),"cm"))
P20 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=Int.P20, fill=SpeciesGroup)) + 
  geom_boxplot() +
  geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "Int.P20"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "Int.P20"])) / 2, linetype="dashed", color = "black") +
  scale_x_discrete(labels = NULL) +
  theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.1,0,0.5),"cm")) +
  ylab("20th percentile of all return intensity") +
  #  ylab(expression(20th~percentile~of~all~return~intensity)) +
  xlab("")
IQ <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=Int.IQ, fill=SpeciesGroup)) + 
  geom_boxplot() +
  geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "Int.IQ"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "Int.IQ"])) / 2, linetype="dashed", color = "black") +
  xlab("Species group") +
  ylab("All return intensity interquartile distance") +
  theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0,0.5,0.1,0.1),"cm"))
PEN <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=First.mean.minus.Last.mean, fill=SpeciesGroup)) + 
  geom_boxplot() +
  geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "First.mean.minus.Last.mean"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "First.mean.minus.Last.mean"])) / 2, linetype="dashed", color = "black") +
  xlab("Species group") +
  ylab("First return mean height - last return mean height (ft)") +
  theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0,0.1,0.1,0.5),"cm"))
grid.arrange(P80, P20, IQ, PEN, nrow = 2, ncol = 2)

# randomly split data into training and test sets
set.seed(seed)
ind = sample(2, nrow(trainData), replace=TRUE, prob=c(0.7,0.3))
trainingData = trainData[ind == 1, ]
testingData = trainData[ind == 2, ]

table(trainingData$SpeciesGroup)
table(testingData$SpeciesGroup)

# *****************************************************************************
# very important!!!
# the row/column arrangement in the confusion matrix output by randomForest
# is the opposite of that output by confusionMatrix
# for the RF output, read across the rows to evaluate the classification
# but for confusionMatrix output, read down the columns to evaluate
# *****************************************************************************

# the code within the "if (TRUE)" or "if (FALSE)" statements is usually only needed
# once. You can either set the statements to "TRUE" and run the code or manually run 
# selected lines using Ctrl-Enter. There are several lines of code that are NOT actually
# needed but were at some point in the model development process. In some cases, lines 
# were commented out but not for all.
#
# set to TRUE to limit the set of variables...makes it much easier to apply the model
# using raster layers
if (FALSE) {
# testing model with fewer variables...logic was to look at variable importance with all variables and 
  # pick the best variable from each "group" of variables. e.g. the full set had Intensity cover above mean, 
  # P70, P75, and P80 in that order so I picked cover above mean and P70.
  #
  # This selection process is subjective but I try to remove metrics that seem redundant. If the variable
  # importance scores give us P70, P75, P90, and P80; it doesn't make sense to include all 4 
  # metrics so I just pick 1 to include. I have done testing (not in this code) and just having 1 of the 4
  # doesn't seem to affect model accuracy.
  #
  # When using the importance scores, correlated predictors cause problems. There is an alternate implementation
  # of RF in the party package that computes conditional permutation importance scores that are
  # supposed to be "better" for use in picking a subset of variables. However, this takes much longer to run
  # so I haven't tried it for this model.
  #
  # I do notice that the random number stream affects the variable importance scores/order. I suppose this
  # means I need to do more testing (enclose the RF call in a testing framework or manually). However, given
  # that I am getting 95%+ prediction accuracy, I'm not sure there is much room for improvement.
  if (modelType == "ALRU") {
    # used these for alder/other model
    trainingData <- trainingData[, c("SpeciesGroup"
                                     , "Int.IQ"
                                     , "Int.P20"
                                     , "First.Int.P80"
                                     , "Int.P30"
                                     , "First.Int.P60"
                                     , "First.mean.minus.Last.mean"
                                     , "First.Int.L.skewness"
                                     , "Last.mean.over.First.mean"
                                     , "Int.P01"
                                     , "X.All.returns.above.4.50.....Total.first.returns....100"
                                     , "Int.RP01"
                                     , "Int.P05"
                                     , "Int.RP05"
                                     , "Int.IQ"
                                     , "Int.L.kurtosis"
    )]
  } else {  
    trainingData <- trainingData[, c("SpeciesGroup"
                                     , "First.Int.RP70"
                                     , "First.Int.L.kurtosis"
                                     , "Percentage.all.returns.above.4.50"
                                     , "First.Canopy.relief.ratio"
                                     , "First.Percentage.all.returns.above.mean"
                                     , "First.Elev.L.skewness"
                                     , "Elev.RP90"
                                     , "First.P90.minus.Last.P10"
    )]
  }
}

# second round of testing with new variable selection idea
# ~94% accuracy 
if (FALSE) {
  set.seed(seed)
  tempData <- trainingData[, c("SpeciesGroup"
                               , "Int.IQ"
                               , "Int.P30"
                               , "First.Int.P60"
                               , "First.mean.minus.Last.mean"
  )]
  # tempData <- trainingData[, c("SpeciesGroup"
  #                                  , "Int.IQ"
  #                                  , "Int.L4"
  #                                  , "Int.L.kurtosis"
  #                                  , "First.Int.skewness"
  #                                  , "Last.mean.over.First.mean"
  # )]
  # tempData <- trainingData[, c("SpeciesGroup"
  #                              , "Int.IQ"
  #                              , "Int.L.kurtosis"
  #                              , "First.Int.P75"
  #                              , "Int.AAD"
  #                              , "First.mean.minus.Last.mean"
  # )]
                               
  tempRF <- randomForest(SpeciesGroup ~ .
                         , data = tempData
                         , importance = (ncol(trainingData) > 10)
                         , sampsize = rep(sum(trainingData$SpeciesGroup == "ALRU")
                                          , nlevels(trainingData$SpeciesGroup)
                         )
                         , mtry = 2, ntree = 1000
  )
  tempRF
  
  typePred <- predict(tempRF, newdata = testingData)
  table(typePred, testingData$SpeciesGroup)
  
  CM <- table(typePred, testingData$SpeciesGroup)
  accuracy <- (sum(diag(CM)))/sum(CM)
  accuracy
  
  m <- confusionMatrix(typePred, testingData$SpeciesGroup)
  m
  
  # do predictions using all data
  typePred <- predict(tempRF, newdata = trainData)
  table(typePred, trainData$SpeciesGroup)
  
  CM <- table(typePred, trainData$SpeciesGroup)
  accuracy <- (sum(diag(CM)))/sum(CM)
  accuracy
  
  m <- confusionMatrix(typePred, trainData$SpeciesGroup)
  m
}

# look at correlation between predictors...not the greatest since we have some highly correlated
# predictors. This is probably why changing mtry to 4 doesn't improve anything
#                               Int.IQ  Int.P20   First.Int.P80 First.mean.minus.Last.mean
#Int.IQ                      1.0000000 -0.5906242     0.8571080                  0.6130366
#Int.P20                    -0.5906242  1.0000000    -0.4417123                 -0.6216737
#First.Int.P80               0.8571080 -0.4417123     1.0000000                  0.5300265
#First.mean.minus.Last.mean  0.6130366 -0.6216737     0.5300265                  1.0000000
#cor(trainingData[, -c(1)])

# use sampsize to balance the sample for each tree. this is needed since we have many more
# conifers than alders and more PSME than TSHE
#
# I tested mtry = 1:4 (and larger values when using all predictors) and mtry=2 seemed to do the best
set.seed(seed)
if (modelType == "ALRU") {
  typeRF <- randomForest(SpeciesGroup ~ .
                       , data = trainingData
                       , importance = (ncol(trainingData) > 10)
                       , sampsize = rep(sum(trainingData$SpeciesGroup == "ALRU")
                          , nlevels(trainingData$SpeciesGroup)
                          )
                       , mtry = 2, ntree = 1000
                       )
} else {
  typeRF <- randomForest(SpeciesGroup ~ .
                         , data = trainingData
                         , importance = TRUE
                         , sampsize = rep(sum(trainingData$SpeciesGroup == "TSHE")
                                          , nlevels(trainingData$SpeciesGroup)
                            )
                         , mtry = 2, ntree = 500
  )
}
typeRF

# alternate way to look at importance...use when we are using the full set of predictors
# relies on importance=TRUE in RF call above
if (ncol(trainingData) > 10) {
  BestVariableCount <- 20
  imp <- data.frame(importance(typeRF, scale=TRUE, type = 2))
  imp$variable <- row.names(imp)
  imp <- imp[order(imp[,1], decreasing = TRUE), ]
  imp <- imp[, c(2, 1)]
  row.names(imp) <- NULL
  v <- imp[c(1:BestVariableCount), 1]
  v
  
  if (FALSE) {
    # look at separability between variables
    for (i in 1:BestVariableCount) {
      t <- data.frame("Variable" = v[i]
                      , "Median ALRU" = median(trainingData[trainingData$SpeciesGroup == "ALRU", v[i]])
                      , "Median CONIFER" = median(trainingData[trainingData$SpeciesGroup == "CONIFER", v[i]])
                      , "RDiff" = abs(median(trainingData[trainingData$SpeciesGroup == "ALRU", v[i]]) - median(trainingData[trainingData$SpeciesGroup == "CONIFER", v[i]])) / max(median(trainingData[, v[i]]))
                      )
      
      if (i == 1)
        results <- t
      else
        results <- rbind(results, t)
    }
    results <- results[order(results[, "RDiff"], decreasing = TRUE), ]
    results
  }
}

# original way to look at variable importance...works with any number of predictors
imp <- as.data.frame(typeRF$importance)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
imp$variable <- row.names(imp)
imp <- imp[, c(5, 3, 4)]
row.names(imp) <- NULL
#imp[c(1:20), ]
imp

# do predictions using the testing data
typePred <- predict(typeRF, newdata = testingData)
table(typePred, testingData$SpeciesGroup)

CM <- table(typePred, testingData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

m <- confusionMatrix(typePred, testingData$SpeciesGroup)
m

# do predictions using all data
typePred <- predict(typeRF, newdata = trainData)
table(typePred, trainData$SpeciesGroup)

CM <- table(typePred, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

m <- confusionMatrix(typePred, trainData$SpeciesGroup)
m

# write off RF model
if (modelType == "ALRU") {
  saveRDS(typeRF, "FINAL_Drone_RF_AlderModel.rds")
} else {
  saveRDS(typeRF, "FINAL_Drone_RF_DW_WH_model.rds")
}

# testing for a simple threshold classifier using Int.IQ (most important predictor with mtry = 1)
# Int.IQ
# result is 92.9% accuracy...not bad but not as good as the RF accuracy
# result is 93.7% accuracy with non-overlapping species
threshold <- (median(trainData$Int.IQ[trainData$SpeciesGroup == "ALRU"]) + median(trainData$Int.IQ[trainData$SpeciesGroup == "CONIFER"])) / 2
RA <- (trainData$Int.IQ >= threshold) * 1
RA[RA == 0] <- 2
RAf <- factor(RA, labels = c("ALRU", "CONIFER"))
table(RAf, trainData$SpeciesGroup)

CM <- table(RAf, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

#Int.P20...alder is less than conifer
# result is 87.3% accuracy
# result is 86.8% accuracy with non-overlapping species
threshold <- (median(trainData$Int.P20[trainData$SpeciesGroup == "ALRU"]) + median(trainData$Int.P20[trainData$SpeciesGroup == "CONIFER"])) / 2
RA <- (trainData$Int.P20 < threshold) * 1
RA[RA == 0] <- 2
RAf <- factor(RA, labels = c("ALRU", "CONIFER"))
table(RAf, trainData$SpeciesGroup)

CM <- table(RAf, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

# First.Int.P80
# result is 86.6% accuracy
# result is 86.5% accuracy with non-overlapping species
threshold <- (median(trainData$First.Int.P80[trainData$SpeciesGroup == "ALRU"]) + median(trainData$First.Int.P80[trainData$SpeciesGroup == "CONIFER"])) / 2
RA <- (trainData$First.Int.P80 >= threshold) * 1
RA[RA == 0] <- 2
RAf <- factor(RA, labels = c("ALRU", "CONIFER"))
table(RAf, trainData$SpeciesGroup)

CM <- table(RAf, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

# First.mean.minus.Last.mean
# result is 88.0% accuracy
# result is 90.8% accuracy with non-overlapping species
threshold <- (median(trainData$First.mean.minus.Last.mean[trainData$SpeciesGroup == "ALRU"]) + median(trainData$Int.IQ[trainData$SpeciesGroup == "CONIFER"])) / 2
RA <- (trainData$First.mean.minus.Last.mean >= threshold) * 1
RA[RA == 0] <- 2
RAf <- factor(RA, labels = c("ALRU", "CONIFER"))
table(RAf, trainData$SpeciesGroup)

CM <- table(RAf, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy


# testing for stability of ntree and mtry parameters
if (TRUE & modelType == "ALRU") {
  for (mtry in seq(1, 16, by = 1)) {
    for (ntree in seq(1000, 5000, by = 1000)) {
      set.seed(seed)
      tRF <- randomForest(SpeciesGroup ~ .
                             , data = trainingData
                             , sampsize = rep(sum(trainingData$SpeciesGroup == "ALRU")
                                              , nlevels(trainingData$SpeciesGroup)
                             )
                             , mtry = mtry, ntree = ntree
      )
      typePred <- predict(tRF, newdata = trainData)
      table(typePred, trainData$SpeciesGroup)
      
      CM <- table(typePred, trainData$SpeciesGroup)
      accuracy <- (sum(diag(CM)))/sum(CM)
      accuracy
      
      t <- data.frame("mtry" = mtry, "ntree" = ntree, "accuracy" = accuracy)
      
      if (mtry == 1 & ntree == 1000)
        result <- t
      else
        result <- rbind(result, t)
    }
  }
  saveRDS(result, "Sensitivity_Drone_RF_AlderModel.rds")

  result <- result[order(result[, "accuracy"], decreasing = TRUE), ]
  result
}
result <-   readRDS("Sensitivity_Drone_RF_AlderModel.rds")












# Code to use R model to map ALRU and PSME over lidar coverage ext --------

# *****************************************************************************
# build a forest mask using cover and P99 height
# we want areas >= 10 feet tall and with 10+% cover
#
# only need to do this once
# *****************************************************************************
if (FALSE) {
  # load layers for mask
  P99 <- raster(paste0(rasterFolder, "elev_P99_4p5plus_5FEET.img"))
  cover <- raster(paste0(rasterFolder, "1st_cover_above4p5_5FEET.img"))
  
  # values for mask...we want >= 10' height and >= 10% cover
  htTreshold <- 10
  coverThreshold <- 10
  
  mask <- (P99 >= htTreshold) & (cover >= coverThreshold)
  writeRaster(mask, paste0(resultsFolder, "ForestMask.img"), format = "HFA", overwrite = TRUE)
  
  #mapview(mask)
  #summary(mask)
  
  rm(P99)
  rm(cover)
}

# *****************************************************************************
# compute various metrics and write to results folder...these are 
# not part of the "regular" metrics...relative to P95
#
# not all of the metrics/layers that are computed or converted to ASCII
# raster format are used in the "final" models. I have code from various 
# testing scenarios that I didn't bother to clean out.
# *****************************************************************************
if (FALSE) {
  IntP95 <- raster(paste0(rasterFolder, "FIRST_RETURNS_int_P95_4p5plus_5FEET.img"))
  IntP01 <- raster(paste0(rasterFolder, "FIRST_RETURNS_int_P01_4p5plus_5FEET.img"))
  IntP05 <- raster(paste0(rasterFolder, "FIRST_RETURNS_int_P05_4p5plus_5FEET.img"))
  IntP70 <- raster(paste0(rasterFolder, "FIRST_RETURNS_int_P70_4p5plus_5FEET.img"))
  
  IntRP01 <- IntP01 / IntP95
  IntRP05 <- IntP05 / IntP95
  IntRP70 <- IntP70 / IntP95
  
  writeRaster(IntRP01, paste0(resultsFolder, "FIRST_RETURNS_int_RP01_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(IntRP01), file= paste0(resultsFolder, "FIRST_RETURNS_int_RP01_4p5plus_5FEET.prj")) 
  writeRaster(IntRP05, paste0(resultsFolder, "FIRST_RETURNS_int_RP05_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(IntRP05), file= paste0(resultsFolder, "FIRST_RETURNS_int_RP05_4p5plus_5FEET.prj")) 
  writeRaster(IntRP70, paste0(resultsFolder, "FIRST_RETURNS_int_RP70_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(IntRP70), file= paste0(resultsFolder, "FIRST_RETURNS_int_RP70_4p5plus_5FEET.prj")) 
  
  rm(IntP01)
  rm(IntP05)
  rm(IntP70)
  rm(IntP95)
  rm(IntRP01)
  rm(IntRP05)
}

# *****************************************************************************
# compute difference between mean of first return heights and mean of last
# return heights
# *****************************************************************************
if (FALSE) {
  firstMean <- raster(paste0(rasterFolder, "FIRST_RETURNS_elev_ave_4p5plus_5FEET.img"))
  lastMean <- raster(paste0(lastRasterFolder, "elev_ave_0p0plus_5FEET.img"))
  firstP90 <- raster(paste0(rasterFolder, "FIRST_RETURNS_elev_P90_4p5plus_5FEET.img"))
  lastP10 <- raster(paste0(lastRasterFolder, "elev_P10_0p0plus_5FEET.img"))
  
  diff <- firstMean - lastMean
  diffP <- firstP90 - lastP10
  
  writeRaster(diff, paste0(resultsFolder, "FirstMean_Minus_LastMean_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(diff), file= paste0(resultsFolder, "FirstMean_Minus_LastMean_4p5plus_5FEET.prj")) 
  writeRaster(diff, paste0(resultsFolder, "FirstMean_Minus_LastMean_4p5plus_5FEET.img"), format = "HFA", overwrite = TRUE, NAflag = -9999)
  writeRaster(diffP, paste0(resultsFolder, "FirstP90_Minus_LastP10_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(diffP), file= paste0(resultsFolder, "FirstP90_Minus_LastP10_4p5plus_5FEET.prj")) 
  writeRaster(diffP, paste0(resultsFolder, "FirstP90_Minus_LastP10_4p5plus_5FEET.img"), format = "HFA", overwrite = TRUE, NAflag = -9999)
  
  rm(firstMean)
  rm(lastMean)
  rm(diff)
  rm(firstP90)
  rm(lastP10)
  rm(diffP)
}

# *****************************************************************************
# layers for AsciiGridPredict need to be in ASCII raster format...do 
# conversion and write to the results folder
# *****************************************************************************
if (FALSE) {
  P90 <- raster(paste0(rasterFolder, "elev_P90_4p5plus_5FEET.img"))
  P95 <- raster(paste0(rasterFolder, "elev_P95_4p5plus_5FEET.img"))
  
  RP90 <- P90 / P95

  writeRaster(RP90, paste0(resultsFolder, "elev_RP90_4p5plus_5FEET.asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
  showWKT(proj4string(RP90), file= paste0(resultsFolder, "elev_RP90_4p5plus_5FEET.prj")) 
  writeRaster(RP90, paste0(resultsFolder, "elev_RP90_4p5plus_5FEET.img"), format = "HFA", overwrite = TRUE, NAflag = -9999)
}

# *****************************************************************************
# layers for AsciiGridPredict need to be in ASCII raster format...do 
# conversion and write to the results folder
# *****************************************************************************
convert <- function(
  layerName,
  rasterFolder,
  resultsFolder
) {
  if (!file.exists(paste0(resultsFolder, layerName, ".asc"))) {
    t <- raster(paste0(rasterFolder, layerName, ".img"))
    writeRaster(t, paste0(resultsFolder, layerName, ".asc"), format = "ascii", overwrite = TRUE, NAflag = -9999)
    showWKT(proj4string(t), file= paste0(resultsFolder, layerName, ".prj")) 
    rm(t)
  } else {
    cat("File already exists\n")
  }
}

if (FALSE) {
  convert("FIRST_RETURNS_int_P60_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("int_P30_4p5plus_5FEET", rasterFolder, resultsFolder)
  
  convert("FIRST_RETURNS_int_P80_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("int_IQ_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("int_P20_4p5plus_5FEET", rasterFolder, resultsFolder)
  
  convert("FIRST_RETURNS_int_kurtosis_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("FIRST_RETURNS_elev_canopy_relief_ratio_5FEET", rasterFolder, resultsFolder)
  convert("all_cover_above4p5_5FEET", rasterFolder, resultsFolder)

  convert("FIRST_RETURNS_all_cover_above_mean_5FEET", rasterFolder, resultsFolder)
  
  # convert("FIRST_RETURNS_all_cover_above_mean_5FEET", rasterFolder, resultsFolder)
  # convert("FIRST_RETURNS_int_P05_4p5plus_5FEET", rasterFolder, resultsFolder)
  # convert("FIRST_RETURNS_int_P01_4p5plus_5FEET", rasterFolder, resultsFolder)
  # convert("FIRST_RETURNS_int_IQ_4p5plus_5FEET", rasterFolder, resultsFolder)
  # convert("FIRST_RETURNS_all_1st_cover_above4p5_5FEET", rasterFolder, resultsFolder)
  # convert("FIRST_RETURNS_int_AAD_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("FIRST_RETURNS_int_Lkurtosis_4p5plus_5FEET", rasterFolder, resultsFolder)
  convert("FIRST_RETURNS_elev_Lskewness_4p5plus_5FEET", rasterFolder, resultsFolder)
}

# *****************************************************************************
# use the model to predict over the area covered by metrics
#
# ***** some of the naming is confusing because the metrics computed for the 
# tree crowns used only first returns. for the cover metrics, we need to also
# use those computed using only first returns. For the tree metrics, the values
# in the FIRST_RETURN layers should match those in the layers without 
# FIRST_RETURN in the names
# *****************************************************************************
if (TRUE) {
  # read RF model
  if (modelType == "ALRU") {
    typeRF <- readRDS("FINAL_RF_AlderModel.rds")
  } else {
    typeRF <- readRDS("FINAL_RF_DW_WH_model.rds")
  }
  
  # build list of input layers...corresponds to the variables used for the classifier...also in the same order
  if (modelType == "ALRU") {
    layers <- list(First.Int.P60 = paste0(resultsFolder, "FIRST_RETURNS_int_P60_4p5plus_5FEET.asc")
                   # , First.Int.P80 = paste0(resultsFolder, "FIRST_RETURNS_int_P80_4p5plus_5FEET.asc")
                   # , Percentage.all.returns.above.mean = paste0(resultsFolder, "FIRST_RETURNS_all_cover_above_mean_5FEET.asc")
                   # , Int.P01 = paste0(resultsFolder, "FIRST_RETURNS_int_P01_4p5plus_5FEET.asc")
                   # , X.All.returns.above.4.50.....Total.first.returns....100 = paste0(resultsFolder, "FIRST_RETURNS_all_1st_cover_above4p5_5FEET.asc")
                   # , Int.RP01 = paste0(resultsFolder, "FIRST_RETURNS_int_RP01_4p5plus_5FEET.asc")
                   # , Int.P05 = paste0(resultsFolder, "FIRST_RETURNS_int_P05_4p5plus_5FEET.asc")
                   # , Int.RP05 = paste0(resultsFolder, "FIRST_RETURNS_int_RP05_4p5plus_5FEET.asc")
                   , Int.IQ = paste0(resultsFolder, "int_IQ_4p5plus_5FEET.asc")
                   # , Int.P20 = paste0(resultsFolder, "int_P20_4p5plus_5FEET.asc")
                   , Int.P30 = paste0(resultsFolder, "int_P30_4p5plus_5FEET.asc")
                   # , Int.L.kurtosis = paste0(resultsFolder, "FIRST_RETURNS_int_Lkurtosis_4p5plus_5FEET.asc")
                   # , Int.AAD = paste0(resultsFolder, "int_AAD_4p5plus_5FEET.asc")
                   , First.mean.minus.Last.mean = paste0(resultsFolder, "FirstMean_Minus_LastMean_4p5plus_5FEET.asc")
    )
    AsciiGridPredict(typeRF, layers, paste0(resultsFolder, "CellType_5FEET.asc"), xtypes = NULL, rows = NULL)
  } else {
    layers <- list(First.Int.RP70 = paste0(resultsFolder, "FIRST_RETURNS_int_RP70_4p5plus_5FEET.asc")
                   , First.Int.L.kurtosis = paste0(resultsFolder, "FIRST_RETURNS_int_Lkurtosis_4p5plus_5FEET.asc")
                   , Percentage.all.returns.above.4.50 = paste0(resultsFolder, "all_cover_above4p5_5FEET.asc")
                   , First.Canopy.relief.ratio = paste0(resultsFolder, "FIRST_RETURNS_elev_canopy_relief_ratio_5FEET.asc")
                   , First.Percentage.all.returns.above.mean = paste0(resultsFolder, "FIRST_RETURNS_all_cover_above_mean_5FEET.asc")
                   , First.Elev.L.skewness = paste0(resultsFolder, "FIRST_RETURNS_elev_Lskewness_4p5plus_5FEET.asc")
                   , Elev.RP90 = paste0(resultsFolder, "elev_RP90_4p5plus_5FEET.asc")
                   , First.P90.minus.Last.P10 = paste0(resultsFolder, "FirstP90_Minus_LastP10_4p5plus_5FEET.asc")
                   # , Int.P01 = paste0(resultsFolder, "FIRST_RETURNS_int_P01_4p5plus_5FEET.asc")
                   # , X.All.returns.above.4.50.....Total.first.returns....100 = paste0(resultsFolder, "FIRST_RETURNS_all_1st_cover_above4p5_5FEET.asc")
                   # , Int.RP01 = paste0(resultsFolder, "FIRST_RETURNS_int_RP01_4p5plus_5FEET.asc")
                   # , Int.P05 = paste0(resultsFolder, "FIRST_RETURNS_int_P05_4p5plus_5FEET.asc")
                   # , Int.RP05 = paste0(resultsFolder, "FIRST_RETURNS_int_RP05_4p5plus_5FEET.asc")
                   # , First.Int.kurtosis = paste0(resultsFolder, "FIRST_RETURNS_int_kurtosis_4p5plus_5FEET.asc")
                   # , Int.AAD = paste0(resultsFolder, "int_AAD_4p5plus_5FEET.asc")
                   # , First.mean.minus.Last.mean = paste0(resultsFolder, "FirstMean_Minus_LastMean_4p5plus_5FEET.asc")
    )
    AsciiGridPredict(typeRF, layers, paste0(resultsFolder, "CellType_DFWH_5FEET.asc"), xtypes = NULL, rows = NULL)
  }
}

# *****************************************************************************
# load the prediction layer and display
#
# the code below mixes the 2 models, applies the mask, and displays things
# *****************************************************************************
# read mask
mask <- raster(paste0(resultsFolder, "ForestMask.img"))

# convert the AsciiGridPredict outputs to IMAGINE format...need to create a projection file
# before doing this or else the IMAGEINE file won't have projection information
if (modelType == "ALRU") {
  type <- raster(paste0(resultsFolder, "CellType_5FEET.asc"))
  showWKT(proj4string(mask), file= paste0(resultsFolder, "CellType_5FEET.prj")) 
  writeRaster(type, paste0(resultsFolder, "CellType_5FEET.img"), format = "HFA", overwrite = TRUE)
} else {
  DFWH_type <- raster(paste0(resultsFolder, "CellType_DFWH_5FEET.asc"))
  showWKT(proj4string(mask), file= paste0(resultsFolder, "CellType_DFWH_5FEET.prj")) 
  writeRaster(DFWH_type, paste0(resultsFolder, "CellType_DFWH_5FEET.img"), format = "HFA", overwrite = TRUE)
}

# reclassify so areas with alder are 0 and areas NOT alder are 1
# we only want results for the DFWH classification for areas that are not alder
if (modelType != "ALRU") {
  nonAlder <- reclassify(type, c(0.5, 1.5, 0,  1.5, 2.5, 1))

  mapview(nonAlder, method = "ngb", na.color = "#FFFFFF80")
  
  DFWH_nonAlder_type <- DFWH_type * nonAlder * mask
  writeRaster(DFWH_nonAlder_type, paste0(resultsFolder, "CellType_DFWH_nonAlder_5FEET.img"), format = "HFA", overwrite = TRUE)
  mapview(DFWH_nonAlder_type, method = "ngb", na.color = "#FFFFFF80")
}

map <- type * mask

#mapview(type, maxpixels =  ncell(type), col.regions = c("green", "yellow"))
mapview(type, zcol = "CellType_5FEET", col.regions = c("green", "yellow"), method = "ngb", na.color = "#FFFFFF80")

# KML layer doesn't look very good due to resampling/compression issues (I think)
type_ll <- projectRaster(type, crs = 4326, method = "ngb")
KML(type_ll, paste0(resultsFolder, "CellType_5FEET.kml"), maxpixels = ncell(type_ll), col = c("green", "yellow"), overwrite = TRUE)
