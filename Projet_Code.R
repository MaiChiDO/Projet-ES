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

#____________________ Importation des données__________________________

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
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, #insig = "blank",
         # Cacher les coefficients de corrélation sur la diagonale
         diag=FALSE)
par(op)

#________________________ Cartographie _________________________________

my.palette <- brewer.pal(n = 7, name = "OrRd")

spplot(gfrance85,"Crime_prop", col.regions = my.palette, cuts = 6, col = "transparent", main="Population par crime\ncontre la propriété")
spplot(gfrance85,"Suicides", col.regions = my.palette, cuts = 6, col = "transparent", main="Population par suicide")
spplot(gfrance85,"Wealth", col.regions = my.palette, cuts = 6, col = "transparent", main="Niveau de richesse")

#_____________________ Matrice de contiguité __________________________

## Definition de voisins
voisins <- poly2nb(gfrance85)

## Calculer des centroides des departements
centroids <- coordinates(gfrance85)

## Matrice de contiguité standardisee
cont.lw <- nb2listw(voisins,style="W")

#____________________ Autocorrélation spatiale _______________________

### Variable expliquée ###

## Test de Moran, Donnees Brutes
moran.test(gfrance85$Suicides, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") # permet de recup la stat de Moran

## Graphique de Moran
op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Suicides,cont.lw,xlab="Population par suicide par département",
           ylab="W * Population par suicide\npar département",
           zero.policy=TRUE)
par(op)

## Test de Geary
geary.test(gfrance85$Suicides, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")


### Variables explicatives ###

## Test de Moran, Donnees Brutes
moran.test(gfrance85$Crime_prop, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") 
moran.test(gfrance85$Wealth, cont.lw, zero.policy=TRUE,randomisation=FALSE,,alternative="two.sided") 

## Graphique de Moran
op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Crime_prop,cont.lw,xlab="Population par crime contre la propriété\npar département",
           ylab="W * Population par crime contre\nla propriété par département",
           zero.policy=TRUE)
par(op)

op <- par(mar=c(5, 6, 4, 2) + 0.1)
moran.plot(x=gfrance85$Wealth,cont.lw,xlab="Niveau de richesse par département",
           ylab="W * Niveau de richesse\npar département",
           zero.policy=TRUE)
par(op)

## Test de Geary
geary.test(gfrance85$Crime_prop, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")
geary.test(gfrance85$Wealth, cont.lw, zero.policy=TRUE,randomisation=FALSE ,alternative="two.sided")

#_________________________ Modèles ___________________________________

## Modèle estimé
modele <- log(Suicides) ~ Crime_prop+Wealth
matrice <- cont.lw

## Modèle MCO
ze.lm <- lm(modele, data=gfrance85)
summary(ze.lm)
## Test de Moran adapté sur les résidus
lm.morantest(ze.lm,matrice,alternative = "two.sided")

## Test LM-Error et LM-Lag
res <- lm.LMtests(ze.lm, matrice, test=c("LMerr", "LMlag",
                                         "RLMerr", "RLMlag"))
summary(res)

## Modèle SEM
ze.sem<-errorsarlm(modele, data=gfrance85, matrice)
summary(ze.sem)
## Test d'Hausman
Hausman.test(ze.sem)

## Modèle SLX
ze.slx <- lmSLX(modele, data=gfrance85, matrice)
summary(ze.slx)

## Modèle SAR
ze.sar<-lagsarlm(modele, data=gfrance85, matrice)
summary(ze.sar)

## Modèle SDM
ze.sardm<-lagsarlm(modele, data=gfrance85, matrice, type="mixed")
summary(ze.sardm)

## Test de l'hypothèse de facteur commun
FC.test<-LR.sarlm(ze.sar,ze.sardm)
print(FC.test)

# _________________________ AIC des modèles ____________________________

aic <- kbl(data.frame("Modèle"="AIC",'SLX'=AIC(ze.slx),
                      'SEM'=AIC(ze.sem),
                      'SAR'=AIC(ze.sar),
                      "SDM"=AIC(ze.sardm)),
           format="latex", align = "cc") %>% kable_classic(full_width = F, html_font = "Cambria")
aic

# _______________ Différents matrices de voisinage _____________________

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

# ____________ SAR pour différents matrices de voisinage _______________

sar.cont <- lagsarlm(modele, data=gfrance85, cont.lw)
sar.2p <- lagsarlm(modele, data=gfrance85, PPV2.lw)
sar.5p <- lagsarlm(modele, data=gfrance85, PPV5.lw)
sar.dist <- lagsarlm(modele, data=gfrance85, dist.lw)

summary(sar.cont)
summary(sar.2p)
summary(sar.5p)
summary(sar.dist)

stargazer(sar.cont,sar.2p,sar.5p,sar.dist,column.labels=c("SAR\nContiguïté","SAR\n2 voisins","SAR\n5 voisins","SAR\nDistance"), align=TRUE, no.space=TRUE)

## Effet marginal
impactSAR <-impacts(sar.cont,listw=cont.lw)
impactSAR



