#!/usr/bin/env Rscript
# Install package if it does not exist yet
if (!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}
require(Cairo)
# This Rscript create a comparision graphic (barplots of sensitivities and accuracies) between the automatic enzyme 
# annotation programs. The Rscript work as an interactive mode as well as in command line. If command line invoked, 
# pass external arguments to R when calling this Rscript which needs, as an input, the comparision files names 
# produce by comparator.py script. The result provides by this program is a file containing the comparison image.

# read file and return a data frame
makeDataFrame <- function(file){
  df = read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE)
  myData = data.frame(prog = df[, c(2)], sp = df[, c(3)])
  return(myData)
}

# return a data frame containing the values of: sensitivity and accuracy of the input given data
get_snp <- function(data) {
  vect <- rep(NA, nrow(data))
  vect <- ifelse(data$prog == 1 & data$sp == 1, "TP", vect)
  vect <- ifelse(data$prog == 1 & data$sp == 0, "FP", vect)
  vect <- ifelse(data$prog == 0 & data$sp == 1, "FN", vect)
  u <- unique(vect)
  l1 = list()
  for (i in 1:length(u)){
    l1[u[i]] = sapply(u[i], function(v) { length(which(vect == v)) })
  }
  sens = as.numeric(l1["TP"]) / (as.numeric(l1["TP"]) + as.numeric(l1["FN"]))
  prec = as.numeric(l1["TP"]) / (as.numeric(l1["TP"]) + as.numeric(l1["FP"]))
  f1score = 2*(sens * prec) / (sens + prec)
  df <- data.frame(sens = sens, prec = prec, f1score = f1score)
  df
}

# build and save the file of comparison between the files passed as arguments
makeBarplot <- function(list_df, output_img, title){
  names = names(list_df)
  Types = c()
  Programs = c()
  Values = c()
  for (i in 1:length(list_df)){
    Types <- c(Types, c("Sensibilité"))
    Programs <- c(Programs, c(names[i]))
    Values <- c(Values, as.numeric(list_df[[i]][1]))
  }
  for (i in 1:length(list_df)){
    Types <- c(Types, c("Précision"))
    Programs <- c(Programs, c(names[i]))
    Values <- c(Values, as.numeric(list_df[[i]][2]))
  }
  for (i in 1:length(list_df)){
    Types <- c(Types, c("F1-Score"))
    Programs <- c(Programs, c(names[i]))
    Values <- c(Values, as.numeric(list_df[[i]][3]))
  }
  df = data.frame(Types = Types, Programs = Programs, Values = Values)
  df$Programs <- as.vector(df$Programs)
  plot <- ggplot(data=df, aes(x=Types, y=Values, fill=Programs)) +
    geom_bar(stat="identity", position=position_dodge(), width = 0.7) +
    xlab("") +
    ylab("") +
    labs(title = paste("Niveau", title)) +
    scale_fill_brewer(palette = "Set2") +
    theme_grey(base_size = 9) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
   #geom_text(aes(label=sprintf("%0.3f", round(Values, digits = 3))), 
  #          position=position_dodge(width=0.25), vjust=0.5, hjust=1.01, size=2.5, angle = 90) +
    theme(
      plot.title   = element_text(color="#000000", size=12),
      legend.title = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(color="#993333", size=8, face="bold.italic", angle = 0),
      axis.title.y = element_text(color="#993333", size=8, face="bold.italic", angle = 90)
    )
 return(plot)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# prédictions du niveau 1
option_list1 = list()
option_list1["PRIAM"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.PRIAM.1.comp"
option_list1["E2P2"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.E2P2.1.comp"
option_list1["KAAS"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.KAAS.1.comp"
option_list1["KOALA"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.KOALA.1.comp"
option_list1["BLAST2GO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.B2G.1.comp"
option_list1["INTERPRO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.IPRSCAN.1.comp"
option_list1["PLEAP_1"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.PLEAP_1.1.comp"
option_list1["PLEAP_0"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev1/tracheophyta.PLEAP_0.1.comp"
progs_names = unlist(names(option_list1))
l1 = list()
for (i in 1:length(progs_names)){
list0 = list()
progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
list0[toupper(progs_names[i])] = list(get_snp(progData))
l1 <- append(l1, list0)
}
p1 = makeBarplot(l1, options$out, "1")
p1 = p1 +theme(legend.position = "none")



# prédictions du niveau 2
option_list1 = list()
option_list1["PRIAM"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.PRIAM.2.comp"
option_list1["E2P2"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.E2P2.2.comp"
option_list1["KAAS"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.KAAS.2.comp"
option_list1["KOALA"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.KOALA.2.comp"
option_list1["BLAST2GO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.B2G.2.comp"
option_list1["INTERPRO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.IPRSCAN.2.comp"
option_list1["PLEAP_1"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.PLEAP_1.2.comp"
option_list1["PLEAP_0"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev2/tracheophyta.PLEAP_0.2.comp"
progs_names = unlist(names(option_list1))
l1 = list()
for (i in 1:length(progs_names)){
list0 = list()
progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
list0[toupper(progs_names[i])] = list(get_snp(progData))
l1 <- append(l1, list0)
}
p2 = makeBarplot(l1, options$out, "2")
p2 = p2 + guides(fill = guide_legend(title = NULL))



# prédictions du niveau 3
option_list1 = list()
option_list1["PRIAM"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.PRIAM.3.comp"
option_list1["E2P2"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.E2P2.3.comp"
option_list1["KAAS"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.KAAS.3.comp"
option_list1["KOALA"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.KOALA.3.comp"
option_list1["BLAST2GO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.B2G.3.comp"
option_list1["INTERPRO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.IPRSCAN.3.comp"
option_list1["PLEAP_1"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.PLEAP_1.3.comp"
option_list1["PLEAP_0"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev3/tracheophyta.PLEAP_0.3.comp"
progs_names = unlist(names(option_list1))
l1 = list()
for (i in 1:length(progs_names)){
list0 = list()
progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
list0[toupper(progs_names[i])] = list(get_snp(progData))
l1 <- append(l1, list0)
}
p3 = makeBarplot(l1, options$out, "3")
p3 = p3 + theme(legend.position = "none")



option_list1 = list()
option_list1["PRIAM"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.PRIAM.4.comp"
option_list1["E2P2"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.E2P2.4.comp"
option_list1["KAAS"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.KAAS.4.comp"
option_list1["KOALA"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.KOALA.4.comp"
option_list1["BLAST2GO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.B2G.4.comp"
option_list1["INTERPRO"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.IPRSCAN.4.comp"
option_list1["PLEAP_1"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.PLEAP_1.4.comp"
option_list1["PLEAP_0"] = "/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/compData_lev4/tracheophyta.PLEAP_0.4.comp"
progs_names = unlist(names(option_list1))
l1 = list()
for (i in 1:length(progs_names)){
list0 = list()
progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
list0[toupper(progs_names[i])] = list(get_snp(progData))
l1 <- append(l1, list0)
}
p4 = makeBarplot(l1, options$out, "4")
p4 = p4 + guides(fill = guide_legend(title = NULL))


##########################
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)
