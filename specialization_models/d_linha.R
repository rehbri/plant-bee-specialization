library(readr)
data.db <- read_delim("bee_d.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

data.dp <- read_delim("plant_d.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)
head(data.db)

hist(data.dp$d)

summary(data.dp$d)



library(car)
library(MuMIn)
library(effects)
library(ggplot2)
library(vegan)
library(lme4) ###lmer & glmer
library(lmerTest)
library(nlme) ##lme & gls
library(MASS) ## glmer.nb & glm.nb
library(rgl)
library(qpcR)
library(corrplot)
library(arm)
library(visreg)
library(DHARMa)
#library(GGdatay)
library(usdm)

#####Model selection?
# Distr de la var de respuesta

require(fitdistrplus)
fit.normal=fitdist(data.db$d,"norm")
par(mfrow=c(1,2))
??fitidist

#plotando as tres distribuições
cdfcomp(fit.normal,horizontals=F, addlegend=T,legendtext=c("Normal"))
qqcomp(fit.normal,addlegend=T,legendtext=c("Normal"))
gofstat(fit.normal)$aic #o menor valor representa a melhor distribuição

###mejor ajuste: normal
par(mfrow=c(1,1))

############ABELHAS
names(data.db)

###USANDO modelos nulos (randomização da abundancia na FD) SEM REMOVER ZEROS e 1s

d_bee<-list()
d_bee[[1]]<-lm(d ~ 1, data=data.db)


d_bee[[2]]<-lm(d ~ ITD+Socialidade+Buzz, data=data.db)
d_bee[[3]]<-lm(d ~ ITD+Socialidade, data=data.db)
d_bee[[4]]<-lm(d ~ ITD+Buzz, data=data.db)
d_bee[[5]]<-lm(d ~ ITD, data=data.db)
d_bee[[6]]<-lm(d ~ Socialidade+Buzz, data=data.db)
d_bee[[7]]<-lm(d ~ Socialidade, data=data.db)
d_bee[[8]]<-lm(d ~ Buzz, data=data.db)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(d_bee, rank="AICc") #Best: Socialidade (1º); Socialidade + Buzz (2º)
summary(d_bee[[7]])
visreg(d_bee[[7]])
summary(d_bee[[6]])






#####

require(fitdistrplus)
fit.normal=fitdist(data.dp$d,"norm")
fit.lognormal=fitdist(data.dp$d,"lnorm")
fit.gamma=fitdist(data.dp$d,"gamma")
par(mfrow=c(1,2))
??fitidist

#plotando as tres distribuições
cdfcomp(fit.normal,horizontals=F, addlegend=T,legendtext=c("Normal"))
qqcomp(fit.normal,addlegend=T,legendtext=c("Normal"))
gofstat(fit.normal)$aic #o menor valor representa a melhor distribuição

###mejor ajuste: normal
par(mfrow=c(1,1))

############PLANTAS
names(data.dp)

###USANDO modelos nulos (randomização da abundancia na FD)

d_plant<-list()
d_plant[[1]]<-lm(d ~ 1, data=data.dp)


d_plant[[2]]<-lm(d ~ Cor+Flower_simetry+Flower_shape+inflorescencia_type, data=data.dp)
d_plant[[3]]<-lm(d ~ Cor+Flower_simetry+Flower_shape, data=data.dp)
d_plant[[4]]<-lm(d ~ Cor+Flower_simetry, data=data.dp)
d_plant[[5]]<-lm(d ~ Cor, data=data.dp)
d_plant[[6]]<-lm(d ~ Flower_simetry+Flower_shape+inflorescencia_type, data=data.dp)
d_plant[[7]]<-lm(d ~ Flower_simetry+Flower_shape, data=data.dp)
d_plant[[8]]<-lm(d ~ Flower_simetry, data=data.dp)
d_plant[[9]]<-lm(d ~ Flower_shape+inflorescencia_type, data=data.dp)
d_plant[[10]]<-lm(d ~ Flower_shape, data=data.dp)
d_plant[[11]]<-lm(d ~ inflorescencia_type, data=data.dp)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(d_plant, rank="AICc") #Nulo

