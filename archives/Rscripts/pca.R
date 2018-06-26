if (!require("ggplot2")) install.packages("ggplot2")

plot_percent_var <- function(pca, pc){
  # Calcule du pourcentage de variance
  percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
  # Préparation d'un tableau avec le numéro des composantes principales 
  # et le pourcentage de variance qui lui est associé
  percent_var_explained <- data.frame(
    PC=1:length(percent_var_explained),
    percent_Var=percent_var_explained
  )
  # Récupérer uniquement le nombre de PC indiqué en argument
  sub_percent_var_explained <- percent_var_explained[1:pc,]
  # Génère le graphique
  p <- ggplot(sub_percent_var_explained, aes(x=PC, y=percent_Var)) + 
    # Génère un barplot
    geom_col()+
    # Utilise le thème "black and white"
    theme_bw() +
    # Renomme l'axe des abscisses
    xlab("PCs") +
    # Renomme l'axe des ordonnées
    ylab("% Variance") +
    # Titre du graphique
    ggtitle("Screeplot")+
    # Option de taille des éléments textuels
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=16),
      legend.text = element_text(size =16),
      legend.title = element_text(size =16 ,face="bold"),
      plot.title = element_text(size=18, face="bold", hjust = 0.5),
      # Astuce pour garder un graphique carré
      aspect.ratio=1
    )
  # Affiche le graphique
  print(p)
}

plot_pca <- function(pca=pca, pc=pc, conditions=conditions, colours=colours){
  # Transforme le nombre de PC en argument en nom de PC 
  PCs <- paste("PC",1:pc, sep="")
  # Calcule le pourcentage de variance par PC
  percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
  # Transforme le vecteur de conditions en un facteur
  cond <- factor(conditions)
  # Crée un autre facteur avec les conditions
  col <- factor(conditions)
  # Change les niveaux du facteur avec la palette de couleur pour attribuer
  # à chaque condition une couleur
  levels(col) <- colours
  # Re-transforme le facteur en vecteur
  col <- as.vector(col)
  # Récupère les scores pour le graphique
  scores <- as.data.frame(pca$x)
  # Génère toutes les combinaisons possibles de PC 
  PCs.combinations <- combn(PCs,2)
  # Génère un graphique pour chaque combinaison
  # avec une boucle apply
  g <- apply(
    PCs.combinations,
    2,
    function(combination)
    {
      p1 <- ggplot(scores, aes_string(x=combination[1], y=combination[2])) +
        # Dessine des points avec une bordure de 0.5 remplis avec une couleur
        geom_point(shape = 21, size = 2.5, stroke=0.5, aes(fill=cond)) +
        # Utilise le thème "black and white"
        theme_bw() +
        # Spécifie la palette de couleur et donne un titre vide à la légende
        scale_fill_manual(
          values=colours,
          name=""
        ) +
        # Renomme le titre des axes des abscisses et des ordonnées en "PCx (pourcentage de variance)" avec 3 chiffres après la virgule
        xlab(paste(combination[1], " (",round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit=3),"%)", sep=""))+
        ylab(paste(combination[2], " (",round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit=3),"%)", sep=""))+
        # Titre du graphique
        ggtitle("PCA")+
        # Option de taille des éléments texte
        theme(
          axis.text=element_text(size=16),
          axis.title=element_text(size=16),
          legend.text = element_text(size =16),
          legend.title = element_text(size =16 ,face="bold"),
          plot.title = element_text(size=18, face="bold", hjust = 0.5),
          # Astuce pour garder un graphique carré
          aspect.ratio=1
        )
      # Affiche le graphique
      print(p1)
    }
  )
}


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

list0 = as.vector(unlist(priamData[1]))
for (i in 1:length(list0)){
  ec0 = unlist(strsplit(list0[i], "\\."))
  if (suppressWarnings(isTRUE(all(ec0 == as.numeric(ec0))))){
    df0 = data.frame(ec = list0[i], prog = "PRIAM", digit1 = as.numeric(ec0[1]), digit2 = as.numeric(ec0[2]), digit3 = as.numeric(ec0[3]), digit4 = as.numeric(ec0[4]))
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

x = df$digit1
y = df$digit2
z = df$digit3
w = df$digit4
data <- data.frame(x, y, z, w)
# Définition des groupes
group <- df$prog
# Définition de la palette de couleur (on peut aussi utiliser RColorBrewer ou tout autre palette déjà faite)
palette <- c("#77b0f3", "#8dcf38", "#fb7072", "#fb7027", "#77b0D3", "#8dff38")

# On lance le calcule de l'ACP
pca <- prcomp(data, center=TRUE, scale=TRUE)
pca
# On affiche le graphique "Screeplot" (pourcentage de variance par composante principale)
plot_percent_var(pca, 4)

# On génère le graphique de l'ACP pour les 2 premières composantes principales
plot_pca(
  pca=pca, 
  pc=2, 
  conditions=group, 
  colours=2:7
)

