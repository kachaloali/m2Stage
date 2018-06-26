#!/usr/bin/env Rscript
# Install package if it does not exist yet
if (!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}

# read file and return a data frame
makeDataFrame <- function(file, prog_name){
  df = read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE)
  prog = rep(prog_name, nrow(df))
  myData = data.frame(ec = df[, c(1)], prog = prog)
  return(myData)
}

priamData = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta.priam_pri.ec.4.list", "PRIAM")
e2p2Data = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta.pf_e2p2.ec.4.list", "E2P2")
kaasData = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta.kaas_kaas.ec.4.list", "KAAS")
koalaData = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta.koala_kaas.ec.4.list", "KOALA")
B2gDatata = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta_b2g.ec.4.list", "Blast2Go")
iprscanDatata = makeDataFrame("/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data.default.param.l4.list/tracheophyta_iprscan.ec.4.list", "IPRSCAN")
df = data.frame(ec = character(), prog = character(), digit1 = numeric(), digit2 = numeric(), digit3 = numeric(), digit4 = numeric())

list0 = as.vector(unlist(priamData[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "PRIAM", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

list0 = as.vector(unlist(e2p2Data[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "E2P2", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

list0 = as.vector(unlist(kaasData[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "KAAS", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

list0 = as.vector(unlist(koalaData[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "KOALA", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

list0 = as.vector(unlist(B2gDatata[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "Blast2Go", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

list0 = as.vector(unlist(iprscanDatata[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "IPRSCAN", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
    df <- rbind(df, df0)
  }
}

program <- df$prog
mesures <- df[,c("digit1", "digit2", "digit3")]
x = df[, c("digit1")]
xlim <- max(x)
y = df[, c("digit2")]
ylim <- max(y)
z = df[, c("digit3")]
zlim <- max(z)
w = df[, c("digit4")]
wlim = max(w)
names(mesures) <- c("digit_1", "digit_2", "digit_3")
#plot(mesures, col = couleur, pch = 19, las = 1)
library(rgl)
plot3d(mesures, type = "s", col = 2:8, size = 2, top = TRUE)
decorate3d(xlim = xlim, ylim = ylim, zlim = zlim, 
           xlab = "x", ylab = "y", zlab = "z", 
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.03)
legend3d("topright", legend = unique(program), pch = 20, col = rainbow(length(unique(program))), cex=.75)















