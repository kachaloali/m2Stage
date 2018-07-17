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
# annotation programs. The Rscript work in command line. If command line invoked, 
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
makeBarplot <- function(list_df, output_img){
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
    scale_fill_brewer(palette = "Set2") +
    theme_grey(base_size = 9) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
    geom_text(aes(label=sprintf("%0.3f", round(Values, digits = 3))), 
              position=position_dodge(width=0.7), vjust=0.5, hjust=1.01, size=2.5, angle = 90) +
    theme(
      plot.title   = element_text(color="#000000", size=12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(color="#993333", size=8, face="bold.italic", angle = 0),
      axis.title.y = element_text(color="#993333", size=8, face="bold.italic", angle = 90)
    )
  # Notice that A4: width=11.69, height=8.27
  width   = 11.69
  height  = 8.27
  ggsave(paste0(output_img, '.png'), plot = plot, path = getwd(), width = width, height = height, dpi = 100, device = "png")
}

# check if the command line is evoked
if (!interactive()){
  # Install package if it does not exist yet
  if (!require("optparse")){
    install.packages("optparse")
    library(optparse)
  }else{
    library(optparse)
  }
  option_list0 = list(
    make_option(c("--progs"), 
                dest = "progs", 
                type = "character", 
                default=NULL, 
                help="The programs files name (Files produce with comparator.py script)", 
                metavar = "character"),
    make_option(c("--out"), 
                dest = "out", 
                type="character", 
                default = "myIMG.png", 
                help="output file name [default= %default]", 
                metavar="character")
  ); 
  opt_parser = OptionParser(option_list=option_list0);
  options = parse_args(opt_parser);
  l0 = unlist(strsplit(options$progs, " ", fixed = TRUE))
  option_list1 = list()
  for (i in seq(from=1, to=length(l0), by=2)){
    option_list1[l0[i]] = l0[i+1]
  }
  progs_names = unlist(names(option_list1))
  l1 = list()
  for (i in 1:length(progs_names)){
    list0 = list()
    progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
    list0[toupper(progs_names[i])] = list(get_snp(progData))
    l1 <- append(l1, list0)
  }
  if ((length(l0) < 2) || (length(l0)%%2 != 0)){
    print_help(opt_parser)
    stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  }else{
    makeBarplot(l1, options$out)
  }
}
