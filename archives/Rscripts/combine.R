library(imager)
library(grid)
library(gridExtra)
library(ggplot2)
fpath1 <- load.image('/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.1.png')
fpath2 <- load.image('/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.2.png')
fpath3 <- load.image('/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.3.png')
fpath4 <- load.image('/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.4.png')



# fpath1 <- '/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.1.png'
# fpath2 <- '/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.2.png'
# fpath3 <- '/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.3.png'
# fpath4 <- '/home/ahassankach/Stage/gitSpace/m2Stage/Rscripts/Data/out075.4.png'

grid.newpage()
grid.draw(ggarrange(fpath1, fpath2, widths = c(2,1)))

library(rvest)
map_il(list(fpath1, fpath2, fpath3, fpath4), load.image) %>% plot




if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(grid)
library(ggpubr)
p1 <- ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) +
  geom_point() 
p2 <- ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) +
  geom_point() + facet_wrap( ~ cyl, ncol=2, scales = "free") +
  guides(colour="none") +
  theme()
grid.newpage()
grid.draw(ggarrange(p1, p2, widths = c(2,1)))