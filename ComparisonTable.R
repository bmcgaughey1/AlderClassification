# code to fit models and generate comparison table for RA classification paper
#
# *****************************************************************************
# packages...may not need all of these for the code below
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
library(kableExtra)
library(coin)
library(extrafont)

# *****************************************************************************
# set up folders using repository structure
# working folder should be the root of the repository
# *****************************************************************************
setwd("~/AlderClassification")

resultsFolder <- "~/AlderClassification/results/"
dataFolder <- "~/AlderClassification/data/"

# *****************************************************************************
# variables to control some things in the code below
# *****************************************************************************
modelType <- "ALRU"
#modeType <- "PSME"   # code not tested gor this model

seed <- 38294579
set.seed(seed)

showInfo <- TRUE

# *****************************************************************************
# read data
# *****************************************************************************
allData <- read.csv(paste0(dataFolder, "AllPlots", "_ModelMetrics.csv"), stringsAsFactors = FALSE)

# make species a factor...necessary for random forest to do classification
allData$Species <- as.factor(allData$Species)

# add new species codes to group conifers and minor species
if (modelType == "ALRU") {
  # create a new "species" code that groups ACCI, PISI, and RHPU into OTHER
  # lump PSME and TSHE into CONIFER type
  allData$SpeciesGroup[allData$Species == "ALRU"] <- "ALRU"
  allData$SpeciesGroup[allData$Species == "PSME"] <- "CONIFER"
  allData$SpeciesGroup[allData$Species == "TSHE"] <- "CONIFER"
  allData$SpeciesGroup[allData$Species == "PISI"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "ACCI"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "RHPU"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "UNK"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "UNKN"] <- "OTHER"
} else {
  # reclass for PSME/TSHE...I tested this and the model has accuracy = ~79%
  # all of the important metrics are height related
  allData$SpeciesGroup[allData$Species == "ALRU"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "PSME"] <- "PSME"
  allData$SpeciesGroup[allData$Species == "TSHE"] <- "TSHE"
  allData$SpeciesGroup[allData$Species == "PISI"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "ACCI"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "RHPU"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "UNK"] <- "OTHER"
  allData$SpeciesGroup[allData$Species == "UNKN"] <- "OTHER"
}

if (showInfo) {
  # show some summary info
  allData %>%
    group_by(Species) %>%
    summarize(
      count = n(),
      aveP70 = mean(Int.P70, na.rm = TRUE),
      aveCovAboveMean = mean(Percentage.all.returns.above.mean, na.rm = TRUE)
    )
}

# drop OTHER group...only has a few individuals with DBH > 10cm
allData <- dplyr::filter(allData, SpeciesGroup != "OTHER")

allData$SpeciesGroup <- as.factor(allData$SpeciesGroup)

# drop trees with anomalies:
# drop trees with values of 1, 2, or 6
# Anomaly_Nu:
# 0: No problem
# 1: Dead tree...drop
# 2: Leaning tree...drop
# 3: Tree shares base with other tree/forked tree
# 4: Tree is on stump
# 5: Tree is less than 0.5 meters away from another tree
# 6: Check in field (location doesn't make sense for some reason or another)...drop
allData <- dplyr::filter(allData, Anomaly_Nu == 0 | Anomaly_Nu == 3 | Anomaly_Nu == 4 | Anomaly_Nu == 5)

# keep trees with DBH >= 10cm
allData <- dplyr::filter(allData, DBH_cm >= 10)

trainData <- allData

# filter using overlap information for crowns
# *********** used for initial testing but found that including the overlapped trees
# gave better predictions for the hold-out testing data
trainData_no_overlap <- dplyr::filter(allData, !Overlapped | !OverlapSpeciesDifferent)

t_absolutely_no_overlap <- dplyr::filter(allData, !Overlapped)

# remove columns for things like return counts
trainData <- trainData[, -c(1:30, 103:104, 118:131, 204:205, 212:217)]
trainData_no_overlap <- trainData_no_overlap[, -c(1:30, 103:104, 118:131, 204:205, 212:217)]

# drop records for trees with too few returns to compute metrics...NA values for most metrics
# this drops 3 trees...1 has only 1 point and the other 2 have 0 points above the height threshold
trainData <- na.omit(trainData)
trainData_no_overlap <- na.omit(trainData_no_overlap)

# look at Mood's median test to check significance of the difference in group median values
if (TRUE) {
  median_test(Int.IQ ~ SpeciesGroup, data = trainData, conf.level = 0.99)
  median_test(First.Int.P60 ~ SpeciesGroup, data = trainData, conf.level = 0.99)
  median_test(First.mean.minus.Last.mean ~ SpeciesGroup, data = trainData, conf.level = 0.99)
  median_test(Int.P30 ~ SpeciesGroup, data = trainData, conf.level = 0.99)
}

# plots for metrics...first uses bse plot() and second uses ggplot2
if (showInfo) {
  # get counts by SpeciesGroup
  table(t_absolutely_no_overlap$SpeciesGroup)
  table(trainData$SpeciesGroup)
  table(trainData_no_overlap$SpeciesGroup)
  
  # show some summary info for variables used for RF model to classify ALRU
  trainData %>%
    group_by(SpeciesGroup) %>%
    summarize(
      Count = n(),
      AveP90 = mean(Elev.P90, na.rm = TRUE),
      FirstIntP80 = mean(First.Int.P80, na.rm = TRUE),
      IntP20 = mean(Int.P20, na.rm = TRUE),
      IntIQ = mean(Int.IQ, na.rm = TRUE),
      FirstLastMeanDiff = mean(First.mean.minus.Last.mean, na.rm = TRUE)
    )
  
  plot(trainData$SpeciesGroup, trainData$Elev.P99, xlab = "Species group", ylab = "P99 height (ft)")
  
  par(mfrow = c(2, 2), oma = c(1, 1, 0, 0), mar = c(2, 4.5, 1, 1))
  plot(trainData$SpeciesGroup, trainData$First.Int.P80, xlab = "Species group", ylab = expression(80^th~percentile~of~first~return~intensity))
  plot(trainData$SpeciesGroup, trainData$Int.P20, xlab = "Species group", ylab = expression(20^th~percentile~of~all~return~intensity))
  plot(trainData$SpeciesGroup, trainData$Int.IQ, xlab = "Species group", ylab = "All return intensity IQ distance")
  plot(trainData$SpeciesGroup, trainData$First.mean.minus.Last.mean, xlab = "Species group", ylab = "First return mean height - last return mean height (ft)")
  par(mfrow = c(1, 1))
  
  # work on improving plots...have problems with alignment due to superscripts in Y-axis labels
  # removed the superscripted parts
  fs <- 20
  ff <- "sans"
  FP60 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=First.Int.P60, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "First.Int.P60"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "First.Int.P60"])) / 2, linetype="dashed", color = "black") +
    #  geom_violin(draw_quantiles = c(.50), trim = FALSE) +
    scale_x_discrete(labels = NULL) +
    ylab("60th percentile of first return intensity") +
    #  ylab(expression(80th~percentile~of~first~return~intensity)) +
    xlab("") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.1,0,0.5),"cm"))
  P30 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=Int.P30, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "Int.P30"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "Int.P30"])) / 2, linetype="dashed", color = "black") +
#    scale_x_discrete(labels = NULL) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.1,0,0.5),"cm")) +
    ylab("30th percentile of all return intensity") +
    #  ylab(expression(20th~percentile~of~all~return~intensity)) +
    xlab("Species group")
  FP80 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=First.Int.P80, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "First.Int.P80"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "First.Int.P80"])) / 2, linetype="dashed", color = "black") +
    #  geom_violin(draw_quantiles = c(.50), trim = FALSE) +
    scale_x_discrete(labels = NULL) +
    ylab("80th percentile of first return intensity") +
    #  ylab(expression(80th~percentile~of~first~return~intensity)) +
    xlab("") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.5,0,0.1),"cm"))
  P20 <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=Int.P20, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "Int.P20"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "Int.P20"])) / 2, linetype="dashed", color = "black") +
    scale_x_discrete(labels = NULL) +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.1,0,0.5),"cm")) +
    ylab("20th percentile of all return intensity") +
    #  ylab(expression(20th~percentile~of~all~return~intensity)) +
    xlab("")
  IQ <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=Int.IQ, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "Int.IQ"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "Int.IQ"])) / 2, linetype="dashed", color = "black") +
    scale_x_discrete(labels = NULL) +
    xlab("") +
    ylab("All return intensity interquartile distance") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.5,0,0.1),"cm"))
  PEN <- ggplot(data = trainData, mapping = aes(x=SpeciesGroup, y=First.mean.minus.Last.mean / 3.2808, fill=SpeciesGroup)) + 
    geom_boxplot() +
    geom_hline(yintercept=(median(trainData[trainData$SpeciesGroup == "ALRU", "First.mean.minus.Last.mean"]) + median(trainData[trainData$SpeciesGroup == "CONIFER", "First.mean.minus.Last.mean"])) / 2 / 3.2808, linetype="dashed", color = "black") +
    xlab("Species group") +
    ylab("Canopy Penetration (m)") +
#    ylab("First return mean height - last return mean height (ft)") +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    theme(legend.position="none", text=element_text(size=fs,  family=ff), plot.margin=unit(c(0.1,0.5,0,0.1),"cm"))
#  grid.arrange(FP80, P20, IQ, PEN, nrow = 2, ncol = 2)
  grid.arrange(IQ, FP60, PEN, P30, nrow = 2, ncol = 2)
  
  # ***** export this plot as a TIFF with width set to 1024 and height to 1572...text size should be OK
  # should be able to automate the save but I had trouble getting things to work...exported manually
}

# function to create confusion matrix and overall accuracy
CM <- function(
  model,
  data
) {
  typePred <- predict(model, newdata = data)
  
  confusionMatrix(typePred, data$SpeciesGroup)
}

# function to compute OOB accuracy for RF model
ACC <- function(
  model
) {
  (100 - model$err.rate[nrow(model$err.rate), 1] * 100)
}

# function to build simple threshold model for discrimination
threshold <- function(
  sourceData,
  applyData,
  variable = "",
  test = "GE"
) {
  threshold <- (median(sourceData[sourceData$SpeciesGroup == "ALRU", variable]) + median(sourceData[sourceData$SpeciesGroup == "CONIFER", variable])) / 2
  if (test == "GE") RA <- (applyData[, variable] >= threshold) * 1
  else RA <- (applyData[, variable] < threshold) * 1
  RA[RA == 0] <- 2
  RAf <- factor(RA, labels = c("ALRU", "CONIFER"))
  table(RAf, applyData$SpeciesGroup)
  
  CM <- table(RAf, applyData$SpeciesGroup)
  accuracy <- (sum(diag(CM)))/sum(CM) * 100
}

seed <- 38294579
testIterations <- 30
set.seed(seed)

# generate seeds for random number generator
seeds <- runif(testIterations) * 2 * 10^9

# fit models testIterations times
for (i in 1:testIterations) {
  seed <- seeds[i]
  cat("Iteration =", i, "-- Seed =", seed, "\n")
  
  # randomly split data into training and test sets
  set.seed(seed)
  ind = sample(2, nrow(trainData), replace=TRUE, prob=c(0.7,0.3))
  trainingData = trainData[ind == 1, ]
  testingData = trainData[ind == 2, ]
  
  set.seed(seed)
  ind = sample(2, nrow(trainData_no_overlap), replace=TRUE, prob=c(0.7,0.3))
  trainingData_no_overlap = trainData_no_overlap[ind == 1, ]
  testingData_no_overlap = trainData_no_overlap[ind == 2, ]
  
  # save data with all metrics for use when evaluating mtry and ntree
  alltrainingData <- trainingData
  alltrainingData_no_overlap <-trainingData_no_overlap
  alltestingData <- testingData
  alltestingData_no_overlap <-testingData_no_overlap
  
  # drop variables...could also do this in the RF call...don't need to do this for the testing data
  trainingData <- trainingData[, c("SpeciesGroup"
                                   , "Int.IQ"
                                   , "Int.P30"
                                   , "First.Int.P60"
                                   , "First.mean.minus.Last.mean"
  )]
  
  trainingData_no_overlap <- trainingData_no_overlap[, c("SpeciesGroup"
                                   , "Int.IQ"
                                   , "Int.P30"
                                   , "First.Int.P60"
                                   , "First.mean.minus.Last.mean"
  )]

  # # drop variables...could also do this in the RF call...don't need to do this for the testing data
  # trainingData <- trainingData[, c("SpeciesGroup"
  #                                  , "Int.IQ"
  #                                  , "Int.P20"
  #                                  , "First.Int.P80"
  #                                  , "First.mean.minus.Last.mean"
  # )]
  # 
  # trainingData_no_overlap <- trainingData_no_overlap[, c("SpeciesGroup"
  #                                                        , "Int.IQ"
  #                                                        , "Int.P20"
  #                                                        , "First.Int.P80"
  #                                                        , "First.mean.minus.Last.mean"
  # )]
  
  # fit RF models
  # RF model with overlapping crowns
  set.seed(seed)
  typeRF <- randomForest(SpeciesGroup ~ .
                         , data = trainingData
                         , importance = FALSE
                         , sampsize = rep(sum(trainingData$SpeciesGroup == "ALRU")
                                          , nlevels(trainingData$SpeciesGroup)
                         )
                         , mtry = 2, ntree = 1000
  )
  typeRF
  
  if (i == 1) {
    results <- data.frame("model" = "RF with overlapping trees"
                          , "Variable" = NA
                          , "IncludesOverlap" = 1
                          , "fitAccuracy" = ACC(typeRF)
                          , "NO_30" = (CM(typeRF, testingData_no_overlap))$overall[1] * 100
                          , "NO_all" = (CM(typeRF, trainData_no_overlap))$overall[1] * 100
                          , "O_30" = (CM(typeRF, testingData))$overall[1] * 100
                          , "O_all" = (CM(typeRF, trainData))$overall[1] * 100
                          , "Kappa" = (CM(typeRF, trainData))$overall[2]
    )
  } else {
    results <- rbind(results, data.frame("model" = "RF with overlapping trees"
                                         , "Variable" = NA
                                         , "IncludesOverlap" = 1
                                         , "fitAccuracy" = ACC(typeRF)
                                         , "NO_30" = (CM(typeRF, testingData_no_overlap))$overall[1] * 100
                                         , "NO_all" = (CM(typeRF, trainData_no_overlap))$overall[1] * 100
                                         , "O_30" = (CM(typeRF, testingData))$overall[1] * 100
                                         , "O_all" = (CM(typeRF, trainData))$overall[1] * 100
                                         , "Kappa" = (CM(typeRF, trainData))$overall[2]
      )
    )
  }
  
  # RF model WITHOUT overlapping crowns
  set.seed(seed)
  typeRF_no_overlap <- randomForest(SpeciesGroup ~ .
                         , data = trainingData_no_overlap
                         , importance = FALSE
                         , sampsize = rep(sum(trainingData_no_overlap$SpeciesGroup == "ALRU")
                                          , nlevels(trainingData_no_overlap$SpeciesGroup)
                         )
                         , mtry = 2, ntree = 1000
  )
  typeRF_no_overlap
  
  results <- rbind(results, data.frame("model" = "RF WITHOUT overlapping trees"
                                       , "Variable" = NA
                                       , "IncludesOverlap" = 0
                                       , "fitAccuracy" = ACC(typeRF_no_overlap)
                                       , "NO_30" = (CM(typeRF_no_overlap, testingData_no_overlap))$overall[1] * 100
                                       , "NO_all" = (CM(typeRF_no_overlap, trainData_no_overlap))$overall[1] * 100
                                       , "O_30" = (CM(typeRF_no_overlap, testingData))$overall[1] * 100
                                       , "O_all" = (CM(typeRF_no_overlap, trainData))$overall[1] * 100
                                       , "Kappa" = (CM(typeRF_no_overlap, trainData))$overall[2]
    )
  )
  
  # RF models with all variables...default mtry, ntree=1000
  set.seed(seed)
  alltypeRF <- randomForest(SpeciesGroup ~ .
                         , data = alltrainingData
                         , importance = FALSE
                         , sampsize = rep(sum(alltrainingData$SpeciesGroup == "ALRU")
                                          , nlevels(alltrainingData$SpeciesGroup)
                         )
                         , ntree = 1000
  )
  alltypeRF
  
  results <- rbind(results, data.frame("model" = "RF all metrics with overlapping trees"
                                       , "Variable" = NA
                                       , "IncludesOverlap" = 1
                                       , "fitAccuracy" = ACC(alltypeRF)
                                       , "NO_30" = (CM(alltypeRF, alltestingData_no_overlap))$overall[1] * 100
                                       , "NO_all" = (CM(alltypeRF, trainData_no_overlap))$overall[1] * 100
                                       , "O_30" = (CM(alltypeRF, alltestingData))$overall[1] * 100
                                       , "O_all" = (CM(alltypeRF, trainData))$overall[1] * 100
                                       , "Kappa" = (CM(alltypeRF, trainData))$overall[2]
  )
  )
  
  set.seed(seed)
  alltypeRF_no_overlap <- randomForest(SpeciesGroup ~ .
                                    , data = alltrainingData_no_overlap
                                    , importance = FALSE
                                    , sampsize = rep(sum(alltrainingData_no_overlap$SpeciesGroup == "ALRU")
                                                     , nlevels(alltrainingData_no_overlap$SpeciesGroup)
                                    )
                                    , ntree = 1000
  )
  alltypeRF_no_overlap
  
  results <- rbind(results, data.frame("model" = "RF all metrics WITHOUT overlapping trees"
                                       , "Variable" = NA
                                       , "IncludesOverlap" = 0
                                       , "fitAccuracy" = ACC(alltypeRF_no_overlap)
                                       , "NO_30" = (CM(alltypeRF_no_overlap, alltestingData_no_overlap))$overall[1] * 100
                                       , "NO_all" = (CM(alltypeRF_no_overlap, trainData_no_overlap))$overall[1] * 100
                                       , "O_30" = (CM(alltypeRF_no_overlap, alltestingData))$overall[1] * 100
                                       , "O_all" = (CM(alltypeRF_no_overlap, trainData))$overall[1] * 100
                                       , "Kappa" = (CM(alltypeRF_no_overlap, trainData))$overall[2]
  )
  )
}

# threshold models...subset of variables
varList <- c("Int.IQ", "Int.P30", "First.Int.P60", "First.mean.minus.Last.mean")

# all variables
varList <- colnames(trainData)[1:203]
for (i in 1:length(varList)) {
  var <- varList[i]

  test <- "GE"
  if (median(trainData[trainData$SpeciesGroup == "ALRU", var]) < median(trainData[trainData$SpeciesGroup == "CONIFER", var]))
    test <- "LT"
  
  # threshold models with overlapping trees
  results <- rbind(results, data.frame("model" = paste("Threshold with overlapping trees:", var)
                                       , "Variable" = var
                                       , "IncludesOverlap" = 1
                                       , "fitAccuracy" = NA
                                       , "NO_30" = NA
                                       , "NO_all" = threshold(trainData, trainData_no_overlap, var, test = test)
                                       , "O_30" = NA
                                       , "O_all" = threshold(trainData, trainData, var, test = test)
                                       , "Kappa" = NA
  )
  )
  
  test <- "GE"
  if (median(trainData_no_overlap[trainData_no_overlap$SpeciesGroup == "ALRU", var]) < median(trainData_no_overlap[trainData_no_overlap$SpeciesGroup == "CONIFER", var]))
    test <- "LT"
  
  # threshold models WITHOUT overlapping trees
  results <- rbind(results, data.frame("model" = paste("Threshold WITHOUT overlapping trees:", var)
                                       , "Variable" = var
                                       , "IncludesOverlap" = 0
                                       , "fitAccuracy" = NA
                                       , "NO_30" = NA
                                       , "NO_all" = threshold(trainData_no_overlap, trainData_no_overlap, var, test = test)
                                       , "O_30" = NA
                                       , "O_all" = threshold(trainData_no_overlap, trainData, var, test = test)
                                       , "Kappa" = NA
  )
  )
}

saveRDS(results, paste0(resultsFolder, "ModelResults_RF_AlderModel.rds"))

results <- readRDS(paste0(resultsFolder, "ModelResults_RF_AlderModel.rds"))

# subset to get all RF models and threshold models for our 4 variables
varList <- c(NA, "Int.IQ", "Int.P30", "First.Int.P60", "First.mean.minus.Last.mean")
t <- results[results$Variable %in% varList, ]

# summarize
sresults <- t %>%
  group_by(model) %>%
  summarize(IncludesOverlap = mean(IncludesOverlap, na.rm = TRUE)
#             , Iterations = n()
            , fitAccuracy = mean(fitAccuracy, na.rm = TRUE)
#            , NO_30 = mean(NO_30, na.rm = TRUE)
            , NO_all = mean(NO_all, na.rm = TRUE)
#            , O_30 = mean(O_30, na.rm = TRUE)
            , O_all = mean(O_all, na.rm = TRUE)
            , Kappa_all = mean(Kappa, na.rm = TRUE)
  )

# manually reorder rows
sresults <- sresults[c(3, 1, 4, 2, 7, 5, 6, 8, 11, 9, 10, 12), ]

# add yes/no to IncludesOverlap column
v <- c("No", "Yes")
sresults$Overlap <- v[sresults$IncludesOverlap + 1]

# add model labels
sresults$ModelName <- c(
  "RF: Variable subset",
  "RF: All variables",
  "RF: Variable subset",
  "RF: All variables",
  "Threshold: Int.IQ",
  "Threshold: First.Int.P60",
  "Threshold: Canopy Penetration",
  "Threshold: Int.P30",
  "Threshold: Int.IQ",
  "Threshold: First.Int.P60",
  "Threshold: Canopy Penetration",
  "Threshold: Int.P30"
)

# rearrange
sresults <- sresults[, c(8, 7, 3, 4, 5, 6)]

# rename columns
colnames(sresults) <- c(
  "Model type",
  "Used overlapped trees",
  "Training data",
  "Non-overlapped trees",
  "All trees",
  "Kappa"
)

# write table
write.csv(sresults, paste0(resultsFolder, "AccuracyTable.csv"), row.names = FALSE, na = "--")
