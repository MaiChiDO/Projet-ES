# Romane LE GOFF, Mai Chi DO, Heloise ROZIER
#_______________________________________________________________________________

rm(list=ls())
#setwd("~/M2 MAS EDP/Econom?trie S & P/ESpatiale/Projet_ES")
#_______________________________________________________________________________

#install.packages("Guerry")
#_______________________________________________________________________________

library(ggplot2)
library(ggmap)
library(sf)
library(leaflet)
library(dplyr)
library(tidyverse)
# Pour les modeles d'econometrie spatiale
library("spdep")
#Pour le calcul de distance (fonction rdist)
library("fields")
#Pour la lecture de table SAS (fonction read.sas7bdat)
library("sas7bdat")
# Pour la cartographie
library("RColorBrewer")
library("classInt")
#_______________________________________________________________________________

library(Guerry)
data(Guerry)
Guerry$Department <- as.character(Guerry$Department)

# Description: 
#Andre-Michel Guerry (1833) collected and analyzed social data on such things as 
#crime,literacy and suicide with the view to determining social laws and the *
#relations among these variables.The Guerry data frame comprises a collection of 
#'moral variables' on the 86 departments of France around 1830.
#A few additional variables have been added from other sources.
#_______________________________________________________________________________

# Cartographie
# Par departement
# Cartes des gros crimes / petits crimes / suicides. Autocorrelation spatiale?


#_________Mise en place de la carte
dpt <- read_sf("dpt")

Guerry$Department <- str_to_upper(Guerry$Department)
# merge seulement les departements qui sont dans guerry : all=FALSE
dp_guerry <-  merge(dpt, Guerry, by.x = "NOM_DEPT", by.y = "Department",all=FALSE)
#76 au lieu de 86, noms de departements qui ont change depuis le temps
nom_guerry <- Guerry$Department %>% sort()
nom_dpt <- dpt$NOM_DEPT %>% sort()

# On travaille en 1833 donc on conserve plutot les anciens noms
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="CHARENTE-MARITIME")] <- "CHARENTE-INFERIEURE" 
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="ALPES-DE-HAUTE-PROVENCE")] <- "BASSES-ALPES" 
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="PYRENEES-ATLANTIQUES")] <- "BASSES-PYRENEES" 
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="SEINE-SAINT-DENIS")] <- "SEINE"
dpt$NOM_DEPT[which(dpt$NOM_DEPT== "VAL-D'OISE" )] <- "SEINE-ET-OISE"
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="SEINE-MARITIME")] <- "SEINE-INFERIEURE"
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="LOIRE-ATLANTIQUE")] <- "LOIRE-INFERIEURE"
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="MEURTHE-ET-MOSELLE")] <- "MEURTHE"
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="COTES-D'ARMOR")] <- "COTES-DU-NORD"
dpt$NOM_DEPT[which(dpt$NOM_DEPT=="TERRITOIRE-DE-BELFORT")] <- "HAUT-RHIN"
#dpt$NOM_DEPT[which(dpt$NOM_DEPT=="HAUTE-CORSE")] <- "CORSE"

dp_guerry <-  merge(dpt, Guerry, by.x = "NOM_DEPT", by.y = "Department",all=F)
# on a maintenant 91 observations. C'est en realite 85 departements. Pb de la seine et seine-et-oise
# Par simplification on pourrait ne pas compter la Corse si jamais on utilise la notion de frontiere

#_________Representation

# gros crimes
ggplot(dp_guerry) + geom_sf(aes(fill=Crime_pers)) +
  scale_fill_continuous(low="sky blue",high="dark blue") +
  ggtitle("Nombre de gros crimes en France (1833)") +
  theme_void()

# petits crimes
ggplot(dp_guerry) + geom_sf(aes(fill=Crime_prop)) +
  scale_fill_continuous(low="sky blue",high="dark blue") +
  ggtitle("Nombre de petits crimes en France (1833)") +
  theme_void()

# suicides
ggplot(dp_guerry) + geom_sf(aes(fill=Suicides)) +
  scale_fill_continuous(low="sky blue",high="dark blue") + 
  ggtitle("Nombre de suicides en France (1833)") +
  theme_void()

# wealth
ggplot(dp_guerry) + geom_sf(aes(fill=Wealth)) +
  scale_fill_continuous(low="sky blue",high="dark blue") + 
  ggtitle("Richesses en France (1833)") +
  theme_void()

# prostitutes => echelle a revoir ?
ggplot(dp_guerry) + geom_sf(aes(fill=Prostitutes)) +
  ggtitle("La prostitution en France (1833)") +
  theme_void()



# Plus grand nombre de suicides dans les departements les plus riches mais moins gros crimes ?
# Il semble en tout cas que l'on va pouvoir mettre des phenomenes en relation
# Lequel ? J'essaye plus loin avec "gros crimes"

# Ou prendre plutot les taux ?
#dp_guerry$tx_crime_pers <- dp_guerry$Crime_pers/sum(dp_guerry$Crime_pers,na.rm=TRUE)

#__________________________________________
### Definition des matrices de voisinage ##
#__________________________________________


### Recuperation des centroides des departements
centro <- st_centroid(dp_guerry$geometry) 
centro <- st_transform(centro,crs=4326)

centro_tab <- do.call(rbind, st_geometry(centro)) %>% 
  as_tibble() %>% 
  setNames(c("lon","lat"))

### Matrice de Contiguite
### Definition des voisins
dpt_sp <- as(dp_guerry, 'Spatial')
voisins <- poly2nb(dpt_sp) #donne l'identifiant de chacun des voisins
summary(voisins)
# L'indiviu 78 est le departement le plus connecte avec 10 liens (Seine et Marne?)
# Les individus 29 et 30 sont la Corse du Sud et du Nord et sont donc voisins entre eux
# Les departements ont en moyenne 5 liens
 
#__________Visualisation des liens

voisins_sf <- as(nb2lines(voisins, coords = coordinates(dpt_sp)), 'sf')
voisins_sf <- st_set_crs(voisins_sf, st_crs(dp_guerry))

ggplot(dp_guerry) + 
  geom_sf(fill = 'slategray1', color = 'white') +
  geom_sf(data = voisins_sf) +
  geom_point(data = centro_tab, aes(x = lon, y = lat)) +
  geom_point(data = centro_tab[78,], aes(x = lon, y = lat),color='red') +
  ggtitle("Visualisation des liens entre departements") +
  coord_sf(crs = 4326) +
  theme_minimal() +
  ylab("Latitude") +
  xlab("Longitude")


#contiguite: frontiere ou pas?
#on passe dune notion de voisinage a une notion de contiguite
### Matrice de contiguite standardisee en ligne, methode par defaut
cont.w <- nb2listw(voisins,style="W")

### Plus proches voisins
### Definition des 2 plus proches voisins
#   k indique le nombre de voisins
cartePPV2.knn <- knearneigh(centro,k=2) 
cartePPV2.nb  <- knn2nb(cartePPV2.knn)
PPV2.w<-nb2listw(cartePPV2.nb,style="W") # Matrice de voisinage des 2 ppv

### Definition des 5 plus proches voisins
cartePPV5.knn<-knearneigh(centro,k=5) 
cartePPV5.nb<- knn2nb(cartePPV5.knn)
PPV5.w<-nb2listw(cartePPV5.nb,style="W")

### Autre option
##basee sur la distance avec des coordonnees ------------------------------------------------------------------------
# wd5 <- dnearneigh(centro, 0, 5)
## en km si true

### Matrice de voisinage fondee sur la distance euclidienne
distance<-rdist(centro_tab,centro_tab) 
# permet de calculer la dist entre tous les individus en fct de leurs coordonnees
# voir aussi fonction "dist" pour d'autres distances

diag(distance)<-0
## rayon de 100 km
distance[distance>=100000 ]<-0
dist<- 1.e6 / (distance*distance) #inverse de la dist au carre
dist[dist>1.e9]<-0
dist.w<-mat2listw(dist, row.names = NULL, style="W")

############################################
### Statistiques descriptives            ###
############################################

### Tableau 1 : Analyses descriptives

library(stargazer)
stargazer(Guerry, type = "text")
# Les variables qui ont des fortes St. Dev. sont interessantes a etudier

############################################
### Graphique de Moran                   ###
############################################

# stat des d'autocorr spatiale

### Test de Moran, Donnees Brutes
moran.test(dp_guerry$Suicides, dist.w, zero.policy=TRUE) # permet de recup la stat de Moran
#I de Moran = ? 0...

### Graphique de Moran
moran.plot(x=dp_guerry$Crime_pers,dist.w,xlab="Criminalite France 1833",
           ylab="Taux dans le voisinage",
           zero.policy=TRUE)




