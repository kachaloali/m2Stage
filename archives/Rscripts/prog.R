library(gplots)
library(gridExtra)
library(ggplot2)
library(grid)
library(pROC)
library(ROCR)
data=read.csv2("/home/ahassankach/Stage/gitSpace/m2Stage/testPriam/Data/e.coli/resCompAnnot.ssp", header=T, sep="")

pred = prediction(data$PriamPred, data[,"spRef"])
perf = performance(pred, "tpr","fpr")
auc <- performance(pred, measure = "auc")@y.values[[1]]
roc.data <- data.frame(specificity=unlist(perf@x.values),
                       sensitivity=unlist(perf@y.values),
                       model="GLM")
basicplot = ggplot(roc.data, aes(x=specificity, ymin=0, ymax=sensitivity)) +
  geom_ribbon(alpha=0.1) +
  geom_line(aes(y=sensitivity))
basicplot + theme(axis.text = element_text(colour = "blue")) + 
  ggtitle("ROC curve") + annotate("text", x = .75, y = .25, label = paste("AUC = ", round(auc, 3)))
auc