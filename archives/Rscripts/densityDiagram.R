# library
if (!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}
if (!require("reshape2")){
  install.packages("reshape2")
  library(reshape2)
}else{
  library(reshape2)
}


blast2go = read.csv("/home/ahassankach/Stage/gitSpace/m2Stage/Data/tracheophyta/tracheophyta_b2g.snp.4.txt", sep = "\t")
kaas = read.csv("/home/ahassankach/Stage/gitSpace/m2Stage/Data/tracheophyta/tracheophyta_kaas.snp.4.txt", sep = "\t")
e2p2 = read.csv("/home/ahassankach/Stage/gitSpace/m2Stage/Data/tracheophyta/tracheophyta.pf_e2p2.snp.4.txt", sep = "\t")
priam = read.csv("/home/ahassankach/Stage/gitSpace/m2Stage/Data/tracheophyta/tracheophyta_pri.snp.4.txt", sep = "\t")
#
par(mfrow=c(1,4))
hist(kaas$sens, breaks = 100, xlim=c(0,1), ylim=c(0, nrow(kaas)), col=rgb(1,0,0,0.5) , xlab="values" , ylab="count" , main="Distribution of sens, spec, prec of kaas program" )
hist(kaas$spec, breaks = 100, xlim=c(0,1) , col=rgb(0,1,0,0.5) , add=T)
hist(kaas$prec, breaks = 100, xlim=c(0,1) , col=rgb(0,0,1,0.5) , add=T)
legend("top", legend=c("sens","spec", "prec"), col=c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
#
hist(blast2go$sens, breaks = 100, xlim=c(0,1), ylim=c(0, nrow(blast2go)), col=rgb(1,0,0,0.5) , xlab="values" , ylab="count" , main="Distribution of sens, spec, prec of blast2go program" )
hist(blast2go$spec, breaks = 100, xlim=c(0,1) , col=rgb(0,1,0,0.5) , add=T)
hist(blast2go$prec, breaks = 100, xlim=c(0,1) , col=rgb(0,0,1,0.5) , add=T)
legend("top", legend=c("sens","spec", "prec"), col=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
#
hist(e2p2$sens, breaks = 100, xlim=c(0,1), ylim=c(0, nrow(e2p2)), col=rgb(1,0,0,0.5) , xlab="values" , ylab="count" , main="Distribution of sens, spec, prec of e2p2 program")
hist(e2p2$spec, breaks = 100, xlim=c(0,1) , col=rgb(0,1,0,0.5) , add=T)
hist(e2p2$prec, breaks = 100, xlim=c(0,1) , col=rgb(0,0,1,0.5) , add=T)
legend("top", legend=c("sens","spec", "prec"), col=c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
#
hist(priam$sens, breaks = 100, xlim=c(0,1), ylim=c(0, nrow(priam)), col=rgb(1,0,0,0.5), xlab="values" , ylab="count" , main="Distribution of sens, spec, prec of priam program")
hist(priam$spec, breaks = 100, xlim=c(0,1) , col=rgb(0,1,0,0.5) , add=T)
hist(priam$prec, breaks = 100, xlim=c(0,1) , col=rgb(0,0,1,0.5) , add=T)
legend("top", legend=c("sens","spec", "prec"), col=c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15 )



# Boxplots
prog = rep("blast2go", nrow(blast2go) + nrow(blast2go) + nrow(blast2go))
prog = c(prog, rep("kaas", nrow(kaas) + nrow(kaas) + nrow(kaas)))
prog = c(prog, rep("e2p2", nrow(e2p2) + nrow(e2p2) + nrow(e2p2)))
prog = c(prog, rep("priam", nrow(priam) + nrow(priam) + nrow(priam)))
prog
Types= c(rep("sensitivity", nrow(blast2go)), rep("specificity", nrow(blast2go)), rep("precision", nrow(blast2go)))
Types= c(Types, rep("sensitivity", nrow(kaas)), rep("specificity", nrow(kaas)), rep("precision", nrow(kaas)))
Types= c(Types, rep("sensitivity", nrow(e2p2)), rep("specificity", nrow(e2p2)), rep("precision", nrow(e2p2)))
Types= c(Types, rep("sensitivity", nrow(priam)), rep("specificity", nrow(priam)), rep("precision", nrow(priam)))
Types
values = c()
values = c(values, blast2go$sens, blast2go$spec, blast2go$prec)
values = c(values, kaas$sens, kaas$spec, kaas$prec)
values = c(values, e2p2$sens, e2p2$spec, e2p2$prec)
values = c(values, priam$sens, priam$spec, priam$prec)
values
data=data.frame(prog, Types , values)
data
# grouped boxplot
ggplot(data, aes(x=prog, y=values, fill=Types, breaks = 100)) + geom_boxplot() +
  labs(fill = "Color type") + 
  theme(
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.title.x = element_text(color="#993333", size=14, face="bold", angle=11),
    axis.title.y = element_text(color="#993333", size=14, face="bold", angle=11)
    
  )
