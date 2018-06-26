# install.packages("magrittr")  # pour le pipe %>%
# install.packages("ade4")      # Calcul de l'ACP
# install.packages("factoextra")# Visualisation  de l'ACP
library(ade4)
library(factoextra)
library(magrittr)
library("factoextra")
library(ggrepel)
# decathlon2.active
# res.pca <- dudi.pca(decathlon2.active,
#                     scannf = FALSE,   # Cacher le scree plot
#                     nf = 5            # Nombre d'axes gardés
# )
# fviz_eig(res.pca)

#####################################################################################################
# res2.pca = dudi.pca(data2,
#                     scannf = FALSE,   # Cacher le scree plot
#                     nf = 5            # Nombre d'axes gardés
# )
# fviz_eig(res2.pca)
#####################################################################################################


# fviz_pca_ind(res2.pca,
#              col.ind = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE
# )


# fviz_pca_var(res2.pca,
#              col.var = "contrib", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     
# )

# 
# fviz_pca_biplot(res2.pca, repel = TRUE,
#                 col.var = "#2E9FDF",
#                 col.ind = "#696969"
# )
data3 = read.table("/home/ahassankach/Stage/gitSpace/m2Stage/scripts/analyzers/ecs_to_pick_l4.tab", header = T, sep = "\t")
column = data3$EC
progs = data3$PROG
names = c()
labels = c()
for (i in 1:nrow(data3)){
  names <- c(names, c(paste(trimws(column[i]), trimws(progs[i]), sep = ".")))
  labels <- c(labels, c(paste(trimws(progs[i]))))
}
row.names(data3) = names
for (j in 2:8){
  data3[,j] = as.numeric(data3[,j])
}
data3 <- data3[-1]

data3 = data3[c("SENS", "PREC")]
data3["labels"] = labels

data3.pca <- prcomp(data3[c("SENS", "PREC")],
                 center = TRUE,
                 scale. = TRUE)
library(ggbiplot)
g <- ggbiplot(data3.pca, 
              obs.scale = 1, 
              var.scale = 1, 
              labels = rownames(data3), 
              labels.size = 1,
              groups = labels, 
              ellipse = TRUE, 
              circle = TRUE)
print(g)


##########################################################
res2.pca = dudi.pca(data3[c("SENS", "PREC")], scannf = FALSE,   # Cacher le scree plot
                    nf = 2           # Nombre d'axes gardés
)
fviz_eig(res2.pca)
fviz_pca_biplot(res2.pca, repel = TRUE,
                col.var = "#2E9FDF",
                col.ind = "#696969"
)
