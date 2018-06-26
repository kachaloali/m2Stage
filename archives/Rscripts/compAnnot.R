# library loading and data preparation: return a data frame object
configAndPrepareData <- function(file){
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(pROC)
  library(ROCR)
  
  setwd(dirname(file))
  dataFrame = read.table(basename(file), header = TRUE, sep = "", stringsAsFactors = FALSE)
  myData = data.frame(progPred = dataFrame[, c(2)], spRef = dataFrame[, c(3)])
  return(myData)
}

# For priam we use the scores instead of binary predictions column
configAndPrepareDataForPriam <- function(file){
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(pROC)
  library(ROCR)
  
  setwd(dirname(file))
  dataFrame = read.table(basename(file), header = TRUE, sep = "", stringsAsFactors = FALSE)
  myData = data.frame(progPred = dataFrame[, c(2)], spRef = dataFrame[, c(3)], scores = dataFrame[, c(4)])
  return(myData)
}


# return a vector containing all of TP, FP, FN and TN
getVector <- function(data) {
  vect <- rep(NA, nrow(data))
  vect <- ifelse(data$progPred == 1 & data$spRef == 1, "TP", vect)
  vect <- ifelse(data$progPred == 1 & data$spRef == 0, "FP", vect)
  vect <- ifelse(data$progPred == 0 & data$spRef == 1, "FN", vect)
  vect <- ifelse(data$progPred == 0 & data$spRef == 0, "TN", vect)
  return(vect)
}

# return a vector containing the proportion of: TP, FP, FN, TN
getCount <- function(vect) {
  u <- unique(vect);
  data.frame(value = u, count = sapply(u, function(v) { length(which(vect == v)) } ))
}

# displays the distribution of: TP, FP, FN, TN
plotDistribution <- function(vect, data){
  ggplot(data=data, aes(x = spRef, y = progPred)) + 
    geom_violin(fill = rgb(1, 1, 1, alpha = 0.6), color = NA) + 
    geom_jitter(aes(color = vect), alpha = 0.6) +
    scale_color_discrete(name = "Types") +
    labs(title=sprintf("Distribution of priam predictions compared to SWISS-PROT annotation using as a standard"))
}

# displays the distribution of: TP, FP, FN, TN in Barplot
printBarplot <- function(vect){
  myData = data.frame(names = vect$value,  value = vect$count)
  ggplot(myData, aes(x = names, y = value)) + geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1")
}


rocdata <- function(dataframe){
  if (length(dataframe[, "response"]) != length(dataframe[, "predictor"])) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(as.factor(dataframe[, "response"]))) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  pred = prediction(dataframe[, c(2)], dataframe[, c(1)])
  perf = performance(pred, "tpr", "fpr")
  roc = data.frame(x = unlist(perf@x.values), y = unlist(perf@y.values))
  roc <- roc[order(roc$x, roc$y),]
  auc <- performance(pred, measure = "auc")
  auc <- auc@y.values[[1]]
  stats <- data.frame (auc = auc)
  return (list(roc = roc, stats = stats))
}

rocplot.multiple <- function(test.data.list, title = "ROC Plot"){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(x)))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc), stats = ldply(plotdata, function(x) x$stats))
  annotation <- with(plotdata$stats, paste("AUC=", signif(auc, 2), sep=""))

  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line() +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = ""))
  return(p)
}


# plot a ROC curve from data frame
plotRocCurves <- function(priamdataframe, dataframe1){
  pred_priam = prediction(priamdataframe$scores, dataframe[, "spRef"])
  perf_priam = performance(pred_priam, "tpr", "fpr")
  auc <- performance(pred_priam, measure = "auc")
  auc <- auc@y.values[[1]]
  roc.data <- data.frame(specificity = unlist(perf_priam@x.values), sensitivity = unlist(perf_priam@y.values), model = "GLM")
  ggplot(roc.data, aes(x = specificity, ymin = 0, ymax = sensitivity)) + geom_ribbon(alpha = 0.1) + 
    geom_line(aes(y = sensitivity)) + ggtitle(paste0("AUC = ", round(auc, 3)))
}



#=================================================================================================================
# prepare data and displays the graphics
priam_dataframe = configAndPrepareDataForPriam("/home/ahassankach/Stage/gitSpace/m2Stage/Data/e.coli/priam_vs_uniprot.txt")
e2p2_dataframe = configAndPrepareData("/home/ahassankach/Stage/gitSpace/m2Stage/Data/e.coli/e2p2_vs_uniprot.txt")

priam = data.frame(response = priam_dataframe$spRef, predictor = priam_dataframe$scores)
e2p2 = data.frame(response = e2p2_dataframe$spRef, predictor = e2p2_dataframe$progPred)
TestData <- list(priam = priam, e2p2 = e2p2, priam1 = priam, priam2 = priam)
plot1 <- rocplot.multiple(TestData, title = "")
plot1
plotdata1 <- llply(TestData1, function(x) with(x, rocdata(x)))
plotdata1



vect = getVector(priamdataframe)
plotDistribution(vect, priamdataframe)

#
vect = getCount(vect)
printBarplot(vect)
barplot(vect$count , density = c(25, 25, 25, 25), angle = c(11, 11, 11, 11) , col = "brown", names.arg = vect$value, axes = TRUE)
axis(2, at = max(vect$count), las = 1)
#
plotRocCurve(dataframe)
