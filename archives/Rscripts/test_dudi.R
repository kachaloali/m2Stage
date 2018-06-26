library(MASS)
data(tdr602)
data(survey)
names(survey)
survey
survey.cc <- survey[complete.cases(survey), ]
survey
sexe <- survey.cc$Sex
couleur <- ifelse(sexe == "Male", rgb(0, 0.5, 1, 0.5), rgb(1, 0.2, 0, 0.5))
couleur
mesures <- survey.cc[,c("Wr.Hnd", "NW.Hnd", "Height")]
names(mesures) <- c("Empan1", "Empan2", "Taille")
plot(mesures, col = couleur, pch = 19, las = 1)

library(rgl)
plot3d(mesures, type = "s", col = couleur)
mesures

library(ade4)
pca1 <- dudi.pca(mesures, scann = FALSE, nf = 3)
names(pca1)
s.label(pca1$li, xax = 1, yax = 2)
gcol = c("blue", "red")
s.class(dfxy = pca1$li, fac = sexe, col = gcol, xax = 1, yax = 2)
s.corcircle(pca1$co, xax = 1, yax = 2)

par(mfrow = c(1, 3))
for (i in 1:3) {
  plot(x = pca1$li[, 1], y = mesures[, i], pch = 19, col = gcol,
       xlab = "First axis of the PCA", las = 1, ylab = colnames(mesures)[i])
}