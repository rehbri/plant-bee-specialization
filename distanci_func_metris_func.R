plant_sp<-read.table("plant_spp.txt",h=T)
bee_sp<-read.table("bee_spp.txt",h=T)



plant_trait<-read.table("plant_trait.txt",h=T, row.names = 1)
bee_trait<-read.table("bee_trait.txt",h=T)
colnames(bee_sp)==rownames(bee_trait)
colnames(plant_sp)==rownames(plant_trait)



library(ade4)
library(vegan)
library(picante)
library(CommEcol)
library(lsr)


#ls()
rm(list=ls())##limpar os documentos existentes
ls()
#removendo MA14 e ST22-08
#DISTANCIA DE GOWER (PAVOINE et al., 2009)

#Distancia Funcional Plantas
names(plant_trait)

#Color:
color<-data.frame(plant_trait[,1])
rownames(color)<-rownames(plant_trait)
colnames(color)<-c("cor")
color

#Simetry
simetry<-data.frame(plant_trait[,2])
rownames(simetry)<-rownames(plant_trait)
colnames(simetry)<-c("simetria")
simetry

#Shape Flower
shape<-data.frame(plant_trait[,3])
rownames(shape)<-rownames(plant_trait)
colnames(shape)<-c("forma")
shape

#inflorrescencia
inflor<-data.frame(plant_trait[,4])
rownames(inflor)<-rownames(plant_trait)
colnames(inflor)<-c("inflorescencia")
inflor


#ARVORE
library(ade4)
ktabplant<-ktab.list.df(list(color, simetry,shape,inflor))

?dist.ktab
func.dist_plant<-dist.ktab(ktabplant,type=c("N","N","N","N"),option=c("scaledBYrange"))
tree_p <- hclust(func.dist_plant, method = "average")#UPGMA agrupamento
plot(tree_p, xlab="", ylab="altura", main="Dendrograma")

library(picante)
tree.plant<-as.phylo(tree_p)
plot(tree.plant)

#### Functional Bees 


names(bee_trait)

#distancia funcional abelhas
#ITD:
itd<-data.frame(bee_trait[,1])
rownames(itd)<-rownames(bee_trait)
colnames(itd)<-c("ITD")
itd

#Socialidade
social<-data.frame(bee_trait[,2])
rownames(social)<-rownames(bee_trait)
colnames(social)<-c("Sociality")
social

#Buzz
buzz<-data.frame(bee_trait[,3:4])
buzz<-prep.binary(buzz, col.blocks = 2, label = "buzz")
rownames(buzz)<-rownames(bee_trait)
colnames(buzz)<-c("Buzz pollination")
buzz



#ARVORE BEES
library(ade4)
ktab_bee<-ktab.list.df(list(itd, social, buzz))


func.dist_bee<-dist.ktab(ktab_bee,type=c("Q","N","B"),option=c("scaledBYrange"))
tree_b <- hclust(func.dist_bee, method = "average")#UPGMA agrupamento
plot(tree_b, xlab="", ylab="altura", main="Dendrograma")

library(picante)
tree.bee<-as.phylo(tree_b)
plot(tree.bee)




#####INDICES DE DIVERSIDADE FUNCIONAL



library(vegan)
library(picante)
library(fBasics)
library(FD)
library(geometry)
library(ade4)



##Distance-Based Functional Diversity Indices (see help of function for more details):
#FRic (Villéger et al. 2008 Ecology):volume
#FEve (Villéger et al. 2008 Ecology):como se distribuem no espaço pesado pela abundância
#FDiv (Villéger et al. 2008 Ecology):
#FDis (Laliberté & Legendre 2010 Ecology):
#Q (Botta-Dukát 2005 JVS):rao, entropia
#FGR (Petchey & Gaston 2006 Ecology Letters): FD do gaston (dendograma)
#CWM (Lavorel et al. 2008 Functional Ecology):

library(FD)
#pode usar a planilha original de traits, mas daí vai usar gower no defaut para calcular a matriz de distÂncia
#pode entrar com a dist calculada com a pavoine, #melhor se tiver variáveis fuzzy e rank#
#defaut usa o melhor número de eixos, pq cada trait é uma coluna
#quando entra com a matriz de distância não calcula o CWM
# só é possível calcular CWM a partir da matriz original de traits

###PLANTAS

DBFD.plantas<-dbFD(func.dist_plant,plant_sp, corr="lingoes")#pode entrar com planilha de abundância se tiver
names(DBFD.plantas)
plant.FD<-cbind.data.frame(DBFD.plantas$nbsp,DBFD.plantas$FRic,DBFD.plantas$FEve,DBFD.plantas$FDis)
colnames(plant.FD)<-c("S","FRic","FEve","FDis")

#Randomizações das amostras: modelos nulos
rand.FD1<-dbFD(func.dist_plant,randomizeMatrix(plant_sp,null.model="frequency"),corr="lingoes")$FRic
rand.FD2<-dbFD(func.dist_plant,randomizeMatrix(plant_sp,null.model="frequency"),corr="lingoes")$FEve
rand.FD3<-dbFD(func.dist_plant,randomizeMatrix(plant_sp,null.model="frequency"), corr="lingoes")$FDis
?randomizeMatrix

plant.null<-cbind.data.frame(rand.FD1,rand.FD2,rand.FD3)
colnames(plant.null)<-c("FRic.null", "FEve.null","FDis.null")
plant.null

FD.plant<-cbind.data.frame(plant.FD,plant.null)
cor(FD.plant[-c(3,4,15,22,24,26,27,28,30),])#removendo os NAs
rownames(FD.plant)


###ABELHAS

DBFD.bees<-dbFD(func.dist_bee,bee_sp, corr="lingoes")#pode entrar com planilha de abundância se tiver
names(DBFD.bees)
bee.FD<-cbind.data.frame(DBFD.bees$nbsp,DBFD.bees$FRic,DBFD.bees$FEve,DBFD.bees$FDis)
colnames(bee.FD)<-c("S","FRic","FEve","FDis")

#Randomizações das amostras: modelos nulos
rand.FD1.bee<-dbFD(func.dist_bee,randomizeMatrix(bee_sp,null.model="independentswap"),corr="lingoes")$FRic
rand.FD2.bee<-dbFD(func.dist_bee,randomizeMatrix(bee_sp,null.model="independentswap"),corr="lingoes")$FEve
rand.FD3.bee<-dbFD(func.dist_bee,randomizeMatrix(bee_sp,null.model="independentswap"), corr="lingoes")$FDis


bee.null<-cbind.data.frame(rand.FD1.bee,rand.FD2.bee,rand.FD3.bee)
colnames(bee.null)<-c("FRic.null", "FEve.null","FDIs.null")
bee.null

FD.bee<-cbind.data.frame(bee.FD,bee.null)
cor(FD.bee[-c(1,17,24,28),])
rownames(FD.bee)








