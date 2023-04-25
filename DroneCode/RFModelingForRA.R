# Red alder classification project with UW using drone lidar data
#

# Introductory comments ---------------------------------------------------
# This is the RF modeling code for the Sappho drone lidar data. This code was
# extracted from a much larger script and simplified (mainly dropped parts that
# weren't needed).

# Setup code -------------------------------------------------------
# *****************************************************************************
# Basic setup stuff
# *****************************************************************************
# I think all of these are used but not all are needed if you don't need 
# fancy graphics or maps
library(sf)
library(raster)
library(terra)
library(mapview)
library(rgdal)
library(rgeos)
library(randomForest)
library(caret)
library(dplyr)
library(ggplot2)
library(gridExtra)

# This is the folder containing the tree data file AllPlots_MOdelMetrics.csv
#
# NOTE: the DroneResults folder is not copied to the GitHub repository
setwd("E:/Backup/R_Stuff/AlderClassification/DroneResults")

# Modeling with RF --------------------------------------------------------

# *****************************************************************************
# fit a classification model
# *****************************************************************************
# this variable controls some of the code below that sets up modeling
# don't use modelType <- "PSME" with the drone lidar data...not sure it will work
modelType <- "ALRU"
#modeType <- "PSME"

seed <- 34567
set.seed(seed)

# read data from above...assumes we are reading the data from the current working folder (set above)
allData <- read.csv("AllPlots_ModelMetrics.csv", stringsAsFactors = FALSE)
trainData <- allData

# make species a factor...necessary for random forest to do classification
trainData$Species <- as.factor(trainData$Species)

# drop trees with anomalies:
# drop dead and leaning tree or trees with code 6
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

# keep trees with DBH >= 10cm
trainData <- dplyr::filter(trainData, DBH_cm >= 10)

# filter using overlap information
# *********** used for initial testing but found that including the overlapped trees
# gave better predictions for the hold-out testing data
#trainData <- dplyr::filter(trainData, !Overlapped | !OverlapSpeciesDifferent)

# get counts by Species
table(trainData$Species)

# lump PSME and TSHE into a CONIFER class
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

# drop OTHER group...only has a few individuals with DBH > 10cm
trainData <- dplyr::filter(trainData, SpeciesGroup != "OTHER")

trainData$SpeciesGroup <- as.factor(trainData$SpeciesGroup)

# get counts by SpeciesGroup...with DBH >= 10cm, we only have 9 OTHER
table(trainData$Species)
table(trainData$SpeciesGroup)

summary(allData$TrueRadius)

# DROP some metrics: mostly return counts and min/max height values
# column 3 has the actual species, column 258 has the species group identifier
trainData <- trainData[, -c(1:30, 103:104, 118:131, 204:205, 212:217)]   # by SpeciesGroup

# drop records for trees with missing metrics
trainData <- na.omit(trainData)

# show some summary info...the focus here is on height. We don't want to see that
# one species is taller than another or else the RF model will tend to use the height
# information over things more related to structure. For example, if PSME is much taller,
# on average, that ALRU, RF will build a model that favors height over most other predictors.
#
# For PSME and ALRU, the heights are about the same so it should be OK to leave the metrics
# closely related to height in the pool of possible metrics.
trainData %>%
  group_by(SpeciesGroup) %>%
  summarize(
    count = n(),
    aveP90 = mean(Elev.P90, na.rm = TRUE),
    aveCovAboveMean = mean(Percentage.all.returns.above.mean, na.rm = TRUE)
  )

# *****************************************************************************
# *****************************************************************************
# at this point, we should have good data


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

# explicitly set the random number seed so we get the same results each time we run the code
set.seed(seed)

# code to limit the variable set used for the classification model and use the model to map species
# over the experimental unit using raster layers with the same metrics
if (TRUE) {
  # subset using just 6 predictors and the species group...this allows us to predict over a larger area without
  # needing all of the metric layers
  trainingDataSubset <- trainingData[, c("Elev.kurtosis", "Elev.skewness", "Int.L2", "First.Elev.kurtosis", "Int.stddev", "Int.L1", "SpeciesGroup")]
  typeRF <- randomForest(SpeciesGroup ~ .
                         , data = trainingDataSubset
                         , importance = (ncol(trainingDataSubset) > 10)
                         , sampsize = rep(sum(trainingDataSubset$SpeciesGroup == "ALRU")
                                          , nlevels(trainingDataSubset$SpeciesGroup)
                         )
                         , mtry = 6, ntree = 1000
  )
  typeRF
  
  # duplicated from above
  library(terra)
  library(raster)
  library(mapview)

  plotList <- c(2, 6, 7, 15, 18, 19, 20, 28, 29, 32, 34, 35)
  for (plot in plotList) {
    # plot <- 34  # plot 8 was not used for model development...presumably no alder
    # resolution <- "5METERS"
    # resolution <- "2p5METERS"    # only have 2.5 raster data for plot 7 & 34
    resolution <- "5METERS"   # only have 1m metrics for 34
    
    # Resolution for the metrics is a bit wonky. Crown circles varied depending on tree height so there isn't a single
    # raster resolution that is "correct". There will be bias if model is built using some resolution and then applied
    # to data at a different resolution but I don't know if this causes problems. Eliminating the bias would be difficult
    # since model was built with different size crown circles.
    
    # read raster layers for metrics and build a raster stack
    # this is a little problematic since the folder names include the date for the AP run and I ran them over 3 days
    # date will be either April 18, 19 or 20 so we can test to see if folder exists for 4/18, if not test 4/19, if not use 4/20
    baseFolder <- paste0("E:/2022_DroneLidar/Sappho_Plot_", plot, "/Products_Sappho_Plot_", plot, "_2023-04-18/FINAL_Sappho_Plot_", plot, "_2023-04-18/Metrics_", resolution, "/")
    if (!file.exists(paste0(baseFolder, paste0("elev_kurtosis_1p37plus_", resolution, ".img")))) {
      baseFolder <- paste0("E:/2022_DroneLidar/Sappho_Plot_", plot, "/Products_Sappho_Plot_", plot, "_2023-04-19/FINAL_Sappho_Plot_", plot, "_2023-04-19/Metrics_", resolution, "/")
      if (!file.exists(paste0(baseFolder, paste0("elev_kurtosis_1p37plus_", resolution, ".img")))) {
        baseFolder <- paste0("E:/2022_DroneLidar/Sappho_Plot_", plot, "/Products_Sappho_Plot_", plot, "_2023-04-20/FINAL_Sappho_Plot_", plot, "_2023-04-20/Metrics_", resolution, "/")
        if (!file.exists(paste0(baseFolder, paste0("elev_kurtosis_1p37plus_", resolution, ".img")))) {
          baseFolder <- paste0("E:/2022_DroneLidar/Sappho_Plot_", plot, "/Products_Sappho_Plot_", plot, "_2023-04-21/FINAL_Sappho_Plot_", plot, "_2023-04-21/Metrics_", resolution, "/")
        }
      }
    }
    
    # build raster stack with layers matching metrics used in RF model
    stack <- c(
      rast(paste0(baseFolder, paste0("elev_kurtosis_1p37plus_", resolution, ".img"))),
      rast(paste0(baseFolder, paste0("elev_skewness_1p37plus_", resolution, ".img"))),
      rast(paste0(baseFolder, paste0("int_L2_1p37plus_", resolution, ".img"))),
      rast(paste0(baseFolder, paste0("FIRST_RETURNS_elev_kurtosis_1p37plus_", resolution, ".img"))),
      rast(paste0(baseFolder, paste0("int_stddev_1p37plus_", resolution, ".img"))),
      rast(paste0(baseFolder, paste0("int_L1_1p37plus_", resolution, ".img")))
    )
    
    t <- rast(paste0(baseFolder, paste0("int_L1_1p37plus_", resolution, ".img")))
    mapview(raster(t))
    
    # assign names that match the names of the lidar metrics used in RF
    names(stack) <- c("Elev.kurtosis", "Elev.skewness", "Int.L2", "First.Elev.kurtosis", "Int.stddev", "Int.L1")
    
    # use our RF model to predict species for the raster stack
    rfp <- predict(stack, typeRF)
    
    # draw a map, class 1 is ALRU, class 2 is CONIFER
    #mapview(raster(rfp), col.regions = c("brown", "green"), at = seq(1, 2, 0.5), alpha.regions = 0.5)
    
    # read trees
    trees <- st_read("AllPlots_Trees.gpkg")
    trees$SpeciesGroup[trees$Species == "ALRU"] <- "ALRU"
    trees$SpeciesGroup[trees$Species == "PSME"] <- "CONIFER"
    trees$SpeciesGroup[trees$Species == "TSHE"] <- "CONIFER"
    trees$SpeciesGroup[trees$Species == "PISI"] <- "OTHER"
    trees$SpeciesGroup[trees$Species == "ACCI"] <- "OTHER"
    trees$SpeciesGroup[trees$Species == "RHPU"] <- "OTHER"
    trees$SpeciesGroup[trees$Species == "UNK"] <- "OTHER"
    trees$SpeciesGroup[trees$Species == "UNKN"] <- "OTHER"
    
    # plot classification results and tree circles.  I can't get this to plot the classification results as categories
    # so I convert and use "special" name so you know the colors and their meaning.
    RAbrown_CONIFERgreen <- raster(rfp)
 
    # clip trees to plot extent
    treesClip <- st_within(trees, sf::st_as_sf(as.polygons(ext(RAbrown_CONIFERgreen), crs = crs(trees))), sparse = FALSE)
    trees <- trees[treesClip, ]

    m1 <- mapview(RAbrown_CONIFERgreen, col.regions = list("brown", "green"), at = seq(1, 2, 0.5), alpha.regions = 0.5, method = "ngb")
    m2 <- mapview(trees, zcol = "SpeciesGroup")
    m <- m1 + m2
  
    # save an image...not rendering the tree
    mapshot(m, file = paste0("Plot_", plot, "_Results.png"))
    
    # write off a KML file for google earth...must reproject to lat-lon
    rfp_ll <- projectRaster(raster(rfp), crs = 4326, method = "ngb")
    KML(rfp_ll, paste0("Plot_", plot, "_CellType_5METERS.kml"), maxpixels = ncell(rfp_ll), col = c("brown", "green"), overwrite = TRUE)
    
    
    
    
    # # load results from Ally's project...in WA state plane south so reprojection is needed (slow)
    # # oldRes <- raster::raster("E:/Backup/R_Stuff/RAClassification/results/CellType_5FEET.img")
    # # oldRes_UTM <- projectRaster(oldRes, crs = 26910, method = "ngb")
    # # oldRes_UTM_crop <- crop(oldRes_UTM, RAbrown_CONIFERgreen)
    # mapview(oldRes_UTM_crop, method = "ngb") +
    #   mapview(RAbrown_CONIFERgreen, maxpixels = 4001 * 4001, col.regions = list("brown", "green"), at = seq(1, 2, 0.5), alpha.regions = 0.5, method = "ngb") +
    #   mapview(trees, zcol = "SpeciesGroup")
    
  }
}

# use sampsize to balance the sample for each tree. this is needed since we have many more
# conifers than alders and more PSME than TSHE
typeRF <- randomForest(SpeciesGroup ~ .
                     , data = trainingData
                     , importance = (ncol(trainingData) > 10)
                     , sampsize = rep(sum(trainingData$SpeciesGroup == "ALRU")
                        , nlevels(trainingData$SpeciesGroup)
                        )
                     , mtry = 6, ntree = 1000
                     )
typeRF

# look at variable importance
imp <- as.data.frame(typeRF$importance)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
imp$variable <- row.names(imp)
imp <- imp[, c(5, 3, 4)]
row.names(imp) <- NULL
imp

# do predictions using the testing data and evaluate accuracy...this is value normally reported
# for model accuracy
typePred <- predict(typeRF, newdata = testingData)
table(typePred, testingData$SpeciesGroup)

CM <- table(typePred, testingData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

m <- confusionMatrix(typePred, testingData$SpeciesGroup)
m

# do predictions using all data...this is not normally done when testing a model
# but it is still useful
typePred <- predict(typeRF, newdata = trainData)
table(typePred, trainData$SpeciesGroup)

CM <- table(typePred, trainData$SpeciesGroup)
accuracy <- (sum(diag(CM)))/sum(CM)
accuracy

m <- confusionMatrix(typePred, trainData$SpeciesGroup)
m

# write off RF model
saveRDS(typeRF, "FINAL_Drone_RF_AlderModel.rds")

# testing for stability of ntree and mtry parameters
# this can be useful to determine the "best" values for ntree and mtry
# this code will take time to run (~10 minutes on my laptop)
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
    typePred <- predict(tRF, newdata = testingData)
    table(typePred, testingData$SpeciesGroup)
    
    CM <- table(typePred, testingData$SpeciesGroup)
    accuracy <- (sum(diag(CM)))/sum(CM)
    accuracy
    
    t <- data.frame("mtry" = mtry, "ntree" = ntree, "accuracy" = accuracy)
    
    if (mtry == 1 & ntree == 1000)
      result <- t
    else
      result <- rbind(result, t)
    
    cat ("mtry = ", mtry, "   ntree = ", ntree, "\n")
  }
}
saveRDS(result, "Sensitivity_Drone_RF_AlderModel.rds")

# Display the results sorted by decreasing accuracy
result <- result[order(result[, "accuracy"], decreasing = TRUE), ]
result
