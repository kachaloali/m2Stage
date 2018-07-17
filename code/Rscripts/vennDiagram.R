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

make_venn_diagram <- function(list_df, level, output_img){
  programs = names(list_df)
  l2 <- list()
  for (i in 1:length(programs)){
    list_ecs = unlist(list_df[programs[i]])
    #print(programs[i])
    #print(length(list_ecs))
    l3 = c()
    for (j in 1:length(list_ecs)){
      #print(list_ecs[j])
      ec_1 = as.character(list_ecs[j])
      ec_1 = strsplit(ec_1, ".", fixed = TRUE)
      ec_1 <- unlist(ec_1)
      remove = c('-')
      ec_1 <- ec_1[! ec_1 %in% remove]
      if (length(ec_1) >= level){
        ec_on_level = ec_1[1:level]
        correct_ec_on_level = paste(ec_on_level, collapse = ".")
        l3 <- c(l3, c(correct_ec_on_level))
      }
    }
    #print(length(l3))
    l2[paste0(toupper(programs[i]))] <- data.frame(l3)
  }
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
    make_option(c("--level"), dest = "level", type = "character", default=NULL, help="The level of ec number classe", metavar = "character"),
    make_option(c("--pleap"), dest = "pleap", type = "character", default=NULL, help="The new pipline file name", metavar = "character"),
    make_option(c("--out"), dest = "out", type="character", default = "myIMG.png", help="output file name [default= %default]", 
                metavar="character")
  ); 
  opt_parser = OptionParser(option_list=option_list);
  options = parse_args(opt_parser);
  if (is.null(options$level)){
    stop("The level of ec number classe is requiered!\n", call. = FALSE)
  }
  l1 = list()
  if (!is.null(options$kaas)){
    kaasData = read.table(options$kaas)
    l1["KAAS"] <- kaasData
  }
  if (!is.null(options$e2p2)){
    e2p2Data = read.table(options$e2p2)
    l1["E2P2"] <- e2p2Data
  }
  if (!is.null(options$priam)){
    priamData = read.table(options$priam)
    l1["PRIAM"] <- priamData
  }
  if (!is.null(options$b2g)){
    b2gData = read.table(options$b2g)
    l1["Blast2go"] <- b2gData
  }
  if (!is.null(options$koala)){
    koalaData = read.table(options$koala)
    l1["KOALA"] <- koalaData
  }
  if (!is.null(options$iprscan)){
    iprscanData = read.table(options$iprscan)
    l1["IPRSCAN"] <- iprscanData
  }
  if (!is.null(options$pleap)){
    pleapData = read.table(options$pleap)
    l1["PLEAP"] <- pleapData
  }
  if (!is.null(options$std)){
    stdData = read.table(options$std)
    l1["STD"] <- stdData
  }
  if (is.null(options$priam) && is.null(options$e2p2) && is.null(options$b2g) && is.null(options$kaas)){
    print_help(opt_parser)
    stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  }else{
   make_venn_diagram(l1, options$level, options$out)
  }
}