if (!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}
library(tidyverse)

# Create the predictor and response variable.
x <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
y <- c(14089,13584,13573,13555,13530,12605,12445,12120,10846,9908,5334)
ylim = 14387
#ylim = ylim + (2000 - ylim%%2000)
data = data.frame(x, y)

# plot graphic
p1 = ggplot(data, aes(x=x, y=y), ylim = ylim) + 
  geom_line(color='#993333') + 
  xlab("F1-Score") +
  ylab("number of annotations") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_y_continuous(limits = c(0, ylim), breaks = seq(0, ylim, by = 2000)) +
  annotate("text", x=0, y=ylim, label= paste0("Total annotation rapported by Swiss-Prot: ", ylim), hjust = 0)
p1



# Create the predictor and response variable.
x <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
y <- c(15053,15053,15053,15053,15053,15053,15053,15045,15011,14970,7067)
ylim = max(y)
lim = 13991
#ylim = ylim + (2000 - ylim%%2000)
data = data.frame(x, y)

# plot graphic
p2 = ggplot(data, aes(x=x, y=y), ylim = ylim) + 
  geom_line(color='#993333') + 
  xlab("F1-Score") +
  ylab("number of annotations") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .1)) +
  scale_y_continuous(limits = c(0, ylim+1000), breaks = seq(0, ylim+1000, by = 2000)) +
  annotate("text", x=0, y=ylim+1000, label= paste0("Total annotation rapported by Swiss-Prot: ", lim), hjust = 0)
p2