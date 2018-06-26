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
    #ggtitle("Performance des prédictions des outils sur les données des trachéophytes") +
    geom_bar(stat="identity", position=position_dodge(), width = 0.25) +
    #xlab("Performance measures") +
    xlab("") +
    ylab("") +
    scale_fill_brewer(palette = "Set2") +
    theme_grey(base_size = 9) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
    geom_text(aes(label=sprintf("%0.3f", round(Values, digits = 3))), 
              position=position_dodge(width=0.25), vjust=0.5, hjust=1.01, size=2.5, angle = 90) +
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
    make_option(c("--progs"), dest = "progs", type = "character", default=NULL, help="Priam file name", metavar = "character"),
    # make_option(c("--priam"), dest = "priam", type = "character", default=NULL, help="Priam file name", metavar = "character"),
    # make_option(c("--e2p2"), dest = "e2p2", type = "character", default=NULL, help="E2P2 file name", metavar = "character"),
    # make_option(c("--b2g"), dest = "b2g", type = "character", default=NULL, help="Blast2go file name", metavar = "character"),
    # make_option(c("--kaas"), dest = "kaas", type = "character", default=NULL, help="KAAS file name", metavar = "character"),
    # make_option(c("--koala"), dest = "koala", type = "character", default=NULL, help="KOALA file name", metavar = "character"),
    # make_option(c("--iprscan"), dest = "iprscan", type = "character", default=NULL, help="InterproScan file name", metavar = "character"),
    # make_option(c("--pleap"), dest = "pleap", type = "character", default=NULL, help="Pleap file name", metavar = "character"),
    make_option(c("--out"), dest = "out", type="character", default = "myIMG.png", help="output file name [default= %default]", 
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
    progData = makeDataFrame(paste0(option_list1[progs_names[i]]))
    progSNP = get_snp(progData)
    print(progs_names[i])
    list0 = list()
    list0[toupper(progs_names[i])] = list(progSNP)
    #print(list0)
    l1 <- append(l1, list0)
  }
  print(l1)
  makeBarplot(l1, options$out)
  stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  l1 = list()
  if (!is.null(options$priam)){
    priamData = makeDataFrame(options$priam)
    priam = get_snp(priamData)
    l1 <- append(l1, list(PRIAM = priam))
  }
  if (!is.null(options$e2p2)){
    e2p2Data = makeDataFrame(options$e2p2)
    e2p2 = get_snp(e2p2Data)
    l1 <- append(l1, list(E2P2 = e2p2))
  }
  if (!is.null(options$b2g)){
    b2gData = makeDataFrame(options$b2g)
    blast2go = get_snp(b2gData)
    l1 <- append(l1, list(Blast2go = blast2go))
  }
  if (!is.null(options$kaas)){
    kaasData = makeDataFrame(options$kaas)
    kaas = get_snp(kaasData)
    l1 <- append(l1, list(KAAS = kaas))
  }
  if (!is.null(options$koala)){
    koalaData = makeDataFrame(options$koala)
    koala = get_snp(koalaData)
    l1 <- append(l1, list(KOALA = koala))
  }
  if (!is.null(options$iprscan)){
    iprscanData = makeDataFrame(options$iprscan)
    iprscan = get_snp(iprscanData)
    l1 <- append(l1, list(IPRSCAN = iprscan))
  }
  if (!is.null(options$pleap)){
    pleapData = makeDataFrame(options$pleap)
    pleap = get_snp(pleapData)
    l1 <- append(l1, list(PLEAP = pleap))
  }
  if (is.null(options$priam) && is.null(options$e2p2) && is.null(options$b2g) && is.null(options$kaas)){
    print_help(opt_parser)
    stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  }else{
    makeBarplot(l1, options$out)
  }
}