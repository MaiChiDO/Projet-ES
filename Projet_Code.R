library(sf)
library(dplyr)
library(tidyverse)
library(spdep)
library(fields)
library(sas7bdat)
library(RColorBrewer)
library(classInt)
library(Guerry)
library(GISTools)
library(sp)
library(stargazer)
library(corrplot)
library(spatialreg)
library(gridExtra)
library(adespatial)
library(spatialreg)
library(kableExtra)

#____________________ Importation des donnees__________________________

data(gfrance85)

#____________________ Statistique descriptive__________________________

df <- data.frame(gfrance85)[,c("Suicides","Crime_prop","Crime_pers","Wealth")]
stargazer(df,type="latex",summary.stat = c("mean", "sd","min","max"))

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
mcor<-cor(df)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(df)
op <- par(mar=c(5, 6, 4, 2) + 0.1)
corrplot(mcor, method="color", col=col(200),
         type="upper", order="hclust",
         addCoef.col = "black", # Ajout du coefficient de correlation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativite
         p.mat = p.mat, sig.level = 0.05, #insig = "blank",
         # Cacher les coefficients de correlation sur la diagonale
         diag=FALSE)
par(op)

#________________________ Cartographie _________________________________

my.palette <- brewer.pal(n = 7, name = "OrRd")

spplot(gfrance85,"Crime_prop", col.regions = my.palette, cuts = 6, col = "transparent", main="Population par crime\ncontre la propriete")
spplot(gfrance85,"Suicides", col.regions = my.palette, cuts = 6, col = "transparent", main="Population par suicide")
spplot(gfrance85,"Wealth", col.regions = my.palette, cuts = 6, col = "transparent", main="Niveau de richesse")

#_____________________ Matrice de contiguite __________________________

## Definition de voisins
voisins <- poly2nb(gfrance85)

## Calculer des centroides des departements
centroids <- coordinates(gfrance85)

## Matrice de contiguite standardisee
cont.lw <- nb2listw(voisins,style="W")

#____________________ Autocorrelation spatiale _______________________

### Variable expliquee ###

## Test de Moran, Donnees Brutes
moran.test(gfrance85$Suicides, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") # permet de recup la stat de Moran

## Graphique de Moran
op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Suicides,cont.lw,xlab="Population par suicide par departement",
           ylab="W * Population par suicide\npar departement",
           zero.policy=TRUE)
par(op)


## Identification des quadrants du diagramme de moran pour chaque observation

#creation d'une nouvelle variable pour identifer les quadrants
gfrance85$quad_sig <- NA
#HH
gfrance85@data[(zx >= 0 & wzx >= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 1.0
#LL
gfrance85@data[(zx <= 0 & wzx <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 2.0
#HL
gfrance85@data[(zx >= 0 & wzx <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 3.0
#LH
gfrance85@data[(zx >= 0 & wzx <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 4.0
#pas significatifs
gfrance85@data[(locm[, 5] > 0.05), "quad_sig"] <- 5.0

#verification
table(gfrance85$quad_sig)

#On affecte 5 aux indicateurs non significatifs
#categorisation
breaks <- seq(1, 5, 1)
# labels des classes 
labels <- c("High-High", "low-Low", "High-Low", "Low-High", "Not Signif.")

# Necessaire pour faire la carte
gfrance85$np <- findInterval(gfrance85$quad_sig, breaks)

# Affectation des couleurs a chaque classe
col.map <- c("darkred", "skyblue2","lightpink","violetred1", "white")
plot(gfrance85, col = col.map[gfrance85$np])  
mtext("Local Moran's I", cex = 1.5, side = 3, line = 1)
legend("topleft", legend = labels, fill = col.map, bty = "n", cex = 0.5)

## Test de Geary
geary.test(gfrance85$Suicides, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")


### Variables explicatives ###

## Test de Moran, Donnees Brutes
moran.test(gfrance85$Crime_prop, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") 
moran.test(gfrance85$Wealth, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") 

## Graphique de Moran
op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Crime_prop,cont.lw,xlab="Population par crime contre la propriete\npar departement",
           ylab="W * Population par crime contre\nla propriete par departement",
           zero.policy=TRUE)
par(op)

op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Wealth,cont.lw,xlab="Niveau de richesse par departement",
           ylab="W * Niveau de richesse\npar departement",
           zero.policy=TRUE)
par(op)

## Test de Geary
geary.test(gfrance85$Crime_prop, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")
geary.test(gfrance85$Wealth, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")

#_________________________ Modeles ___________________________________

## Modele estime
modele <- Suicides ~ Crime_prop + Wealth
matrice <- cont.lw

## Modele MCO
ze.lm <- lm(modele, data=gfrance85)
summary(ze.lm)
## Test de Moran adapte sur les residus
lm.morantest(ze.lm,matrice,alternative = "two.sided")

## Test LM-Error et LM-Lag
res <- lm.LMtests(ze.lm, matrice, test=c("LMerr", "LMlag",
                                         "RLMerr", "RLMlag"))
summary(res)

## Modele SEM
ze.sem<-errorsarlm(modele, data=gfrance85, matrice)
summary(ze.sem)
## Test d'Hausman
Hausman.test(ze.sem)

## Modele SLX
ze.slx <- lmSLX(modele, data=gfrance85, matrice)
summary(ze.slx)

## Modele SAR
ze.sar<-lagsarlm(modele, data=gfrance85, matrice)
summary(ze.sar)

## Modele SDM
ze.sardm<-lagsarlm(modele, data=gfrance85, matrice, type="mixed")
summary(ze.sardm)

## Test de l'hypothese de facteur commun
FC.test<-LR.sarlm(ze.sar,ze.sardm)
print(FC.test)

stargazer(ze.lm,ze.sar,ze.sardm,ze.sem,ze.slx,column.labels=c("MCO","SAR","SDM","SEM","SLX"), align=TRUE, no.space=TRUE)

# _________________________ AIC des modèles ____________________________

aic <- kbl(data.frame("Modele"="AIC",'SLX'=AIC(ze.slx),
                      'SEM'=AIC(ze.sem),
                      'SAR'=AIC(ze.sar),
                      "SDM"=AIC(ze.sardm)),
           format="latex", align = "cc") %>% kable_classic(full_width = F, html_font = "Cambria")
aic

# _______________ Differents matrices de voisinage _____________________

## Definition des 2 plus proches voisins
##   k indique le nombre de voisins
cartePPV2.knn <- knearneigh(centroids,k=2) 
cartePPV2.nb  <- knn2nb(cartePPV2.knn)
PPV2.lw<-nb2listw(cartePPV2.nb,style="W") # Matrice de voisinage des 2 ppv

## Definition des 5 plus proches voisins
cartePPV5.knn<-knearneigh(centroids,k=5) 
cartePPV5.nb<- knn2nb(cartePPV5.knn)
PPV5.lw<-nb2listw(cartePPV5.nb,style="W")

## Matrice de voisinage fondee sur la distance euclidienne
distance<-rdist(centroids) 
# permet de calculer la dist entre tous les individus en fct de leurs coordonnees
# voir aussi fonction "dist" pour d'autres distances

diag(distance)<-0
## rayon de 100 km
distance[distance>=100000 ]<-0
dist<- 1.e6 / (distance*distance) #inverse de la dist au carre
dist[dist>1.e9]<-0
dist.lw<-mat2listw(dist, row.names = NULL, style="W")

# ____________ SAR pour differents matrices de voisinage _______________

sar.cont <- lagsarlm(modele, data=gfrance85, cont.lw)
sar.2p <- lagsarlm(modele, data=gfrance85, PPV2.lw)
sar.5p <- lagsarlm(modele, data=gfrance85, PPV5.lw)
sar.dist <- lagsarlm(modele, data=gfrance85, dist.lw)

summary(sar.cont)
summary(sar.2p)
summary(sar.5p)
summary(sar.dist)

stargazer(sar.cont,sar.2p,sar.5p,sar.dist,column.labels=c("SAR\nContiguete","SAR\n2 voisins","SAR\n5 voisins","SAR\nDistance"), align=TRUE, no.space=TRUE)

## Effet marginal
impactSAR <-impacts(sar.cont, listw=cont.lw)
impactSAR
