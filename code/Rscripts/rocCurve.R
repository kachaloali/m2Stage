calculate_roc <- function(df, title){
  df$std <- as.factor(df$std)
  df$pred <- as.factor(df$pred)
  
  tpr <- function(df, threshold){
    sum(df$score >= threshold & df$std == 1) / sum(df$std == 1)
  }
  
  fpr <- function(df, threshold) {
    sum(df$score >= threshold & df$std == 0) / sum(df$std == 0)
  }
  roc <- data.frame(threshold = seq(0,1,length.out=length(df$score)), tpr = NA, fpr = NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  return (ggplot(roc, aes(fpr, tpr)) + 
    geom_line(color = "#993333") +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "longdash", color = "#993333") + 
    labs(title = sprintf(title)) + 
    xlab("1-Précision") + 
    ylab("Sensibilité") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .125)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .125)) +
    annotate("text", x = 0.30, y = 0.55, label= sprintf("auc: %0.2f", auc(df$std, df$score)), 
             color = "#993333", hjust = 0, angle = 0)
    )
}


# f1 = "/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/tracheophyta.new.matrix.1.txt"
# f2 = "/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/tracheophyta.new.matrix.2.txt"
# f3 = "/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/tracheophyta.new.matrix.3.txt"
# f4 = "/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/tracheophyta.new.matrix.4.txt"
# 
# 
# df1 <- read.table(f1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# df2 <- read.table(f2, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# df3 <- read.table(f3, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# df4 <- read.table(f4, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # roc1 <- calculate_roc(df1, "Niveau 1")
# # roc2 <- calculate_roc(df2, "Niveau 2")
# roc3 <- calculate_roc(df3, "Niveau 3")
# roc4 <- calculate_roc(df4, "Niveau 4")
# library(gridExtra)
# #grid.arrange(roc1, roc2, roc3, roc4, ncol = 2, nrow = 2) # if wanting to plot 4 curves at the same time
# grid.arrange(roc3, roc4)


# check if the command line is evoked
if (!interactive()){
  # Install package if it does not exist yet
  if (!require("optparse")){
    install.packages("optparse")
    library(optparse)
  }else{
    library(optparse)
  }
  library(ggplot2)
  library(pROC)
  library(gridExtra)
  option_list0 = list(
    make_option(c("--f"), 
                dest = "files", 
                type = "character", 
                default=NULL, 
                help="The roc curves files name (Files produce with make_matrix.py script)", 
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
  l0 = unlist(strsplit(options$files, " ", fixed = TRUE))
  
  roclist = list()
  for (i in 1:length(l0)){
    df <- read.table(l0[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    roc1 = calculate_roc(df, "")
    roclist[[i]] = roc1
  }

  if (length(l0) < 1){
    print_help(opt_parser)
    stop("At least one argument must be supplied, Please try again\n", call. = FALSE)
  }else{
    n <- length(l0)
    nCol <- floor(sqrt(n))
    print(nCol)
    plot = grid.arrange(unlist(roclist), ncol=nCol)
    ggsave(paste0(options$out, '.png'), plot = plot, path = getwd(), width = width, height = height, dpi = 100, device = "png")
    print("finished")
  }
}
