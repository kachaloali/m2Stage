#!/usr/bin/env Rscript
# Install package if it does not exist yet
if (!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}
if (!require(futile.logger)){
  install.packages("futile.logger")
  library(futile.logger)
}else{
  library(futile.logger)
}
flog.threshold(ERROR)
# making venn diagram
library(VennDiagram)

# This Rscript compute a venn diagram graphic from a list of ec predicted by the automated annotation programs 
# The Rscript work as an interactive mode. If command line invoked, pass external arguments to R when calling 
# this Rscript which needs, as an input, the list of ec files names produce by get_ec_list.py script. 
# The result provides by this program is a file containing the venn diagram image.
# read file and return a data frame
makeDataFrame <- function(file){
  df = read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE)
  myData = data.frame(ec = df[, c(1)], pred = df[, c(2)], std = df[, c(3)])
  return(myData)
}

make_venn_diagram <- function(list_df, ecNum, output_img){
  programs = names(list_df)
  l2 <- list()
  ec_to_list = unlist(strsplit(ecNum, ".", fixed = TRUE))
  level = length(ec_to_list)
  #print(typeof(ecNum))
  for (i in 1:length(programs)){
    #list_ecs = unlist(list_df[programs[i]])
    print (programs[i])
    data = list_df[programs[i]]
    df <- data.frame(data, byrow = T, stringsAsFactors=FALSE)
    #vect <- ifelse(data$pred == 1 & data$std == 1, "TP", vect)
    #print(length(vect))
    #print(nrow(df))
    l3 = c()
    std_list = c()
    k = 0
    l = 0
    m = 0
    n = 0
    for (j in 1:nrow(df)){
      pred = as.character(df[j,][[2]])
      std =  as.character(df[j,][[3]])
      ec = unlist(strsplit(as.character(df[j,][[1]]), ".", fixed = TRUE))
      ec <- ec[! ec %in% c('-')]
      if (length(ec) >= level){
        ec_on_level = ec[1:level]
        correct_ec_on_level = paste(ec_on_level, collapse = ".")
        if (correct_ec_on_level == ecNum){
          if (pred == "1" & std == "1"){
            l3 <- c(l3, c(paste0("TP", k)))
            k = k + 1
          }else if (pred == "1" & std == "0"){
            l3 <- c(l3, c(paste0(programs[i], l)))
            l = l + 1
          }else if (pred == "0" & std == "1"){
            l3 <- c(l3, c(paste0(programs[i], l)))
            l = l + 1
          }
          
          if (i == 1){
            if (std == "1"){
              std_list <- c(std_list, c(paste0("TP", m)))
              m = m + 1
            }else if (std == "0"){
              std_list <- c(std_list, c(paste0("STD", n)))
              n = n + 1
            }
          }
        }
      }
    }
    if (i == 1){
      l2[paste0(toupper("std"))] <- data.frame(std_list)
    }
    if(length(l3) > 0){
      l2[paste0(toupper(programs[i]))] <- data.frame(l3)
    }
  }
  if (length(l2) > 0){
    plot <- venn.diagram(l2, 
                         fill = 1:length(l2), alpha = 0.3, filename = NULL,
                         output = FALSE , compression = "lzw", category = names(l2),
                         cat.fontface = "bold", cat.fontfamily = "time", lwd = "2", 
                         main = "Diagramme de Venn de prédiction des numéros ECs", main.fontface = "bold"
    )
    # Notice that A4: width=11.69, height=8.27
    width   = 11.69
    height  = 8.27
    ggsave(paste0(output_img, '.png'), plot = plot, path = getwd(), width = width, height = height, dpi = 100, device = "png")
  }else{
    stop("The ec number passed is not detected by any program\n", call. = FALSE)
  }
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
  option_list = list(
    make_option(c("--priam"), dest = "priam", type = "character", default=NULL, help="Priam file name", metavar = "character"),
    make_option(c("--e2p2"), dest = "e2p2", type = "character", default=NULL, help="E2P2 file name", metavar = "character"),
    make_option(c("--b2g"), dest = "b2g", type = "character", default=NULL, help="Blast2go file name", metavar = "character"),
    make_option(c("--kaas"), dest = "kaas", type = "character", default=NULL, help="KAAS file name", metavar = "character"),
    make_option(c("--koala"), dest = "koala", type = "character", default=NULL, help="KOALA file name", metavar = "character"),
    make_option(c("--iprscan"), dest = "iprscan", type = "character", default=NULL, help="InterproScan file name", metavar = "character"),
    make_option(c("--std"), dest = "std", type = "character", default=NULL, help="The swiss prot standard file name", metavar = "character"),
    make_option(c("--ec"), dest = "ec_number", type = "character", default=NULL, help="The level of ec number classe", metavar = "character"),
    make_option(c("--pleap"), dest = "pleap", type = "character", default=NULL, help="The new pipline file name", metavar = "character"),
    make_option(c("--out"), dest = "out", type="character", default = "myIMG.png", help="output file name [default= %default]", 
                metavar="character")
  ); 
  opt_parser = OptionParser(option_list=option_list);
  options = parse_args(opt_parser);
  if (is.null(options$ec_number)){
    stop("The ec number is requiered!\n", call. = FALSE)
  }
  l1 = list()
  if (!is.null(options$kaas)){
    kaasData = makeDataFrame(options$kaas)
    l1 <- append(l1, list(KAAS = kaasData))
    #kaasData = read.table(options$kaas)
    #l1["KAAS"] <- kaasData
  }
  if (!is.null(options$e2p2)){
    e2p2Data = makeDataFrame(options$e2p2)
    l1 <- append(l1, list(E2P2 = e2p2Data))
    #e2p2Data = read.table(options$e2p2)
    #l1["E2P2"] <- e2p2Data
  }
  if (!is.null(options$priam)){
    priamData = makeDataFrame(options$priam)
    l1 <- append(l1, list(PRIAM = priamData))
    #priamData = read.table(options$priam)
    #l1["PRIAM"] <- priamData
  }
  if (!is.null(options$b2g)){
    b2gData = makeDataFrame(options$b2g)
    l1 <- append(l1, list(Blast2go = b2gData))
    #b2gData = read.table(options$b2g)
    #l1["Blast2go"] <- b2gData
  }
  if (!is.null(options$koala)){
    koalaData = makeDataFrame(options$koala)
    l1 <- append(l1, list(KOALA = koalaData))
    #koalaData = read.table(options$koala)
    #l1["KOALA"] <- koalaData
  }
  if (!is.null(options$iprscan)){
    iprscanData = makeDataFrame(options$iprscan)
    l1 <- append(l1, list(IPRSCAN = iprscanData))
    #iprscanData = read.table(options$iprscan)
    #l1["IPRSCAN"] <- iprscanData
  }
  if (!is.null(options$pleap)){
    pleapData = makeDataFrame(options$pleap)
    l1 <- append(l1, list(PLEAP = pleapData))
    #pleapData = read.table(options$pleap)
    #l1["PLEAP"] <- pleapData
  }
  # if (!is.null(options$std)){
  #   stdData = makeDataFrame(options$std)
  #   #stdData = read.table(options$std)
  #   l1["STD"] <- stdData
  # }
  if (is.null(options$priam) && is.null(options$e2p2) && is.null(options$b2g) && is.null(options$kaas)){
    print_help(opt_parser)
    stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  }else{
   make_venn_diagram(l1, options$ec_number, options$out)
  }
}