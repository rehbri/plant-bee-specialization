
library(readr)
data <- read_delim("data_FD_ID_IE_novo.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

head(data)

hist(data$ID)

summary(data$ID)

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
library(GGdatay)
library(usdm)

#####Model selection?
# Distr de la var de respuesta

require(fitdistrplus)
fit.normal=fitdist(data$ID,"norm")
fit.lognormal=fitdist(data$ID,"lnorm")
fit.gamma=fitdist(data$ID,"gamma")
par(mfrow=c(1,2))

#plotando as tres distribuições
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição

###mejor ajuste: normal
par(mfrow=c(1,1))

#

ID_Bee<-list()
ID_Bee[[1]]<-lm(ID ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
ID_Bee[[2]]<-lm(ID ~ FDis.ab+multi_riq, data=data)
ID_Bee[[3]]<-lm(ID ~ FRic.ab+multi_riq,data=data)
ID_Bee[[4]]<-lm(ID ~ FEve.ab+multi_riq, data=data)

ID_Bee[[5]]<-lm(ID ~ FDiv.ab+multi_riq, data=data)
ID_Bee[[6]]<-lm(ID ~ FRic.pa+multi_riq, data=data)
ID_Bee[[7]]<-lm(ID ~ FDis.pa+multi_riq, data=data)
ID_Bee[[8]]<-lm(ID ~ FEve.pa+multi_riq, data=data)
ID_Bee[[9]]<-lm(ID ~ FDiv.pa+multi_riq, data=data)

ID_Bee[[10]]<-lm(ID ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_Bee, rank="AICc") #Best: NULL

 head(model.sel(ID_Bee, rank="AICc"),4)


summary(ID_Bee[[4]])

visreg(ID_Bee[[8]],"FDis.pa")

visreg(ID_Bee[[8]],by="multi_riq","FDis.pa",overlay=T)

visreg(ID_Bee[[8]],by="multi_riq","FDis.pa",overlay=T,breaks=c(24,480))

visreg(ID_Bee[[8]],cond=list(multi_riq=100),"FDis.pa",overlay=T)#definir valor para o tamanho da rede

####plot fdiv.pa

visreg(ID_Bee[[10]],"FDiv.pa")

visreg(ID_Bee[[10]],by="multi_riq","FDiv.pa",overlay=T)

visreg(ID_Bee[[10]],by="multi_riq","FDiv.pa",overlay=T,breaks=c(24,480))

visreg(ID_Bee[[10]],cond=list(multi_riq=100),"FDiv.pa",overlay=T)#definir valor para o tamanho da rede


###SÓ ABUNDANCIA

ID_Bee<-list()
ID_Bee[[1]]<-lm(ID ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
ID_Bee[[2]]<-lm(ID ~ FDis.ab+multi_riq, data=data)
#ID_Bee[[3]]<-lm(ID ~ FRic.ab+multi_riq,data=data)
ID_Bee[[3]]<-lm(ID ~ FEve.ab+multi_riq, data=data)
ID_Bee[[4]]<-lm(ID ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_Bee, rank="AICc") #Best: NULL
summary(ID_Bee[[1]])



##SÓ PRESENÇA E SUSENCIA
ID_Bee<-list()
ID_Bee[[1]]<-lm(ID ~ 1, data=data)
#ID_Bee[[2]]<-lm(ID ~ FRic.pa+multi_riq, data=data)
ID_Bee[[2]]<-lm(ID ~ FDis.pa+multi_riq, data=data)
ID_Bee[[3]]<-lm(ID ~ FEve.pa+multi_riq, data=data)
ID_Bee[[4]]<-lm(ID ~ multi_riq, data=data)


model.sel(ID_Bee, rank="AICc") #Best: NULL



summary(ID_Bee[[1]])

visreg(ID_Bee[[2]],"FDis.pa")

visreg(ID_Bee[[2]],by="multi_riq","FDis.pa",overlay=T)

visreg(ID_Bee[[2]],by="multi_riq","FDis.pa",overlay=T,breaks=c(24,480))

visreg(ID_Bee[[2]],cond=list(multi_riq=100),"FDis.pa",overlay=T,xlab="Dispersão funcional das abelhas (FDis)",ylab="Diversidade de interações (ID)")#definir valor para o tamanho da rede

visreg(ID_Bee[[4]],cond=list(multi_riq=100),"multi_riq",overlay=T,xlab="Tamanho da rede (Nsize)", ylab="Diversidade de interações (ID)")#definir valor para o tamanho da rede



library(visreg) # Certifique-se de que você já instalou e carregou a biblioteca visreg

?visreg

summary(data$multi_riq)


####para IE_sum

#####Model selection?
# Distr de la var de respuesta
hist(data$IE_sum)
require(fitdistrplus)
fit.normal=fitdist(data$IE_sum,"norm")
fit.lognormal=fitdist(data$IE_sum,"lnorm")
fit.gamma=fitdist(data$IE_sum,"gamma")
par(mfrow=c(1,2))

#plotando as tres distribuições
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição
#SE EU QUISER USAR O LOG OU O GAMA?
summary(data$IE_sum)
###mejor ajuste: normal
par(mfrow=c(1,1))

#aBUNDANCIA
IE_sum_Bee<-list()
IE_sum_Bee[[1]]<-lm(IE_sum ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
IE_sum_Bee[[2]]<-lm(IE_sum ~ FDis.ab+multi_riq, data=data)
IE_sum_Bee[[3]]<-lm(IE_sum ~ FEve.ab+multi_riq, data=data)
IE_sum_Bee[[4]]<-lm(IE_sum ~ multi_riq, data=data)
model.sel(IE_sum_Bee, rank="AICc") #Best: NULL
summary(IE_sum_Bee[[3]])

#####PRESENÇA E AUSENCIA

IE_sum_Bee<-list()
IE_sum_Bee[[1]]<-lm(IE_sum ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
IE_sum_Bee[[2]]<-lm(IE_sum ~ FDis.pa+multi_riq, data=data)
IE_sum_Bee[[3]]<-lm(IE_sum ~ FEve.pa+multi_riq, data=data)
IE_sum_Bee[[4]]<-lm(IE_sum ~ multi_riq, data=data)
model.sel(IE_sum_Bee, rank="AICc") #Best: NULL

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(IE_sum_Bee, rank="AICc") #Best: NULL

summary(IE_sum_Bee[[2]])

library(car)
car::Anova(IE_sum_Bee[[4]])
