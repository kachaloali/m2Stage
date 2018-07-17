data1 = read.table("/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/out.heatmap.sens.3.txt", header = T, sep = "\t", dec = ",")
data2 =  read.table("/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/out.heatmap.prec.3.txt", header = T, sep = "\t", dec = ",")
data3 = read.table("/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/out.heatmap.f1score.3.txt", header = T, sep = "\t", dec = ",")

### Retrive the ec numbers at the first culumn
column1 = data1$EC
column2 = data2$EC
column3 = data3$EC

### Will be used as names to the matrix taking by the heatmap function
names1 = c()
names2 = c()
names3 = c()
###
for (i in 1:nrow(data1)){
  names1 <- c(names1, c(trimws(column1[i])))
}
for (i in 1:nrow(data2)){
  names2 <- c(names2, c(trimws(column2[i])))
}
for (i in 1:nrow(data3)){
  names3 <- c(names3, c(trimws(column3[i])))
}


### Making ECs numbers as the names of the matrix
row.names(data1) = names1
row.names(data2) = names2
row.names(data3) = names3
data1


### Convert to numeric, the different values of performances mesures
for (j in 2:7){
  data1[,j] = as.numeric(data1[,j])
}
for (j in 2:7){
  data2[,j] = as.numeric(data2[,j])
}
for (j in 2:7){
  data3[,j] = as.numeric(data3[,j])
}

### To remove the first culmn, that's to say, the culumn of the ECs numbers.
data1 <- data1[-1]
data2 <- data2[-1]
data3 <- data3[-1]


### Convert all data as matrix
data1 = as.matrix(data1)
data2 = as.matrix(data2)
data3 = as.matrix(data3)


### Names of the programs used to analyze the data set of reference
x1 = colnames(data1)
y1 = rownames(data1)
x2 = colnames(data2)
y2 = rownames(data2)
x3 = colnames(data3)
y3 = rownames(data3)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(RColorBrewer) # Rcolorbrewer palette: is interesting for setting color of the heatmap image
library(ComplexHeatmap) # Is a library that allows it to build a heatmap
# heat1 = Heatmap(data1[, 1:10], col = colorRampPalette(brewer.pal(8, "Reds"))(10), name = "heat1")
# heat2 = Heatmap(data2[, 1:10], col = colorRampPalette(brewer.pal(8, "Reds"))(10), name = "heat2")
# heat3 = Heatmap(data3, col = colorRampPalette(brewer.pal(8, "Reds"))(10), name = "heat3")
# heat_list = heat1 + heat2
# draw(heat_list, row_title = "Two heatmaps, row title", row_title_gp = gpar(col = "red"),
#      column_title = "Two heatmaps, column title", column_title_side = "bottom")


library(d3heatmap) # Is a library that allows it to build a heatmap too.
# d3heat1 = d3heatmap(data1, symm=FALSE, scale = "none", col = colorRampPalette(brewer.pal(8, "Reds"))(10))
# d3heat2 = d3heatmap(data2, symm=FALSE, scale = "none", col = colorRampPalette(brewer.pal(8, "Reds"))(10))
# d3heat3 = d3heatmap(data3, symm=FALSE, scale = "none", col = colorRampPalette(brewer.pal(8, "Reds"))(10))



library(iheatmapr) # Is a library that allows it to build a heatmap
# for plotting one heatmap
# main_heatmap = main_heatmap(data1, name = "Sens", colors =colorRampPalette(brewer.pal(8, "Reds"))(10)) %>%
#   add_col_labels(ticktext = x1, font = list(size = 10)) %>%
#   add_row_labels(size = 0.3,font = list(size = 6)) %>%
#   add_row_clustering() %>% 
#   add_col_clustering()
#   main_heatmap


#%%%%%%%%%%%%%%%%
library(iheatmapr)
# for plotting multiple heatmap
main_heatmap = main_heatmap(data1, name = "Sens", colors =colorRampPalette(brewer.pal(8, "Reds"))(10)) %>%
  add_col_labels(ticktext = x1, font = list(size = 10)) %>%
  add_row_labels(size = 0.3,font = list(size = 6)) %>%
  add_row_clustering() %>% 
  add_col_clustering() %>%
add_main_heatmap(data2, name = "Prec", colors =colorRampPalette(brewer.pal(8, "Blues"))(10)) %>%
  add_col_labels(ticktext = x2, font = list(size = 10)) %>%
  add_col_clustering() %>%
add_main_heatmap(data3, name = "F1-Score", colors =colorRampPalette(brewer.pal(8, "Greens"))(10))  %>%
  add_col_labels(ticktext = x3, font = list(size = 10)) %>%
  add_col_clustering()
main_heatmap
