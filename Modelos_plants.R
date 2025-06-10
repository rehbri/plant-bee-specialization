library(readr)
setwd("C:\\Users\\renat\\OneDrive - Universidade Federal do Pará - UFPA\\MESTRADO\\Projeto\\dados\\script_correção_ID\\projeto_correcao_ID\\FD")
data <- read_delim("data_FD_ID_IE_PLANTS.csv",  delim = ";", escape_double = FALSE, trim_ws = TRUE)

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

ID_Plant<-list()
ID_Plant[[1]]<-lm(ID ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
ID_Plant[[2]]<-lm(ID ~ FDis.ab+multi_riq, data=data)
ID_Plant[[3]]<-lm(ID ~ FRic.ab+multi_riq,data=data)
ID_Plant[[4]]<-lm(ID ~ FEve.ab+multi_riq, data=data)
#ID_Plant[[5]]<-lm(ID ~ RaoQ+multi_riq, data=data)
ID_Plant[[5]]<-lm(ID ~ FDiv.ab+multi_riq, data=data)
ID_Plant[[6]]<-lm(ID ~ FRic.pa+multi_riq, data=data)
ID_Plant[[7]]<-lm(ID ~ FDis.pa+multi_riq, data=data)
ID_Plant[[6]]<-lm(ID ~ FEve.pa+multi_riq, data=data)
ID_Plant[[7]]<-lm(ID ~ FDiv.pa+multi_riq, data=data)
ID_Plant[[8]]<-lm(ID ~ multi_riq, data=data)


###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_Plant, rank="AICc") #Best: NULL

head(model.sel(ID_Plant, rank="AICc"),5)


summary(ID_Plant[[3]])

visreg(ID_Plant[[5]],"RaoQ")

visreg(ID_Plant[[5]],by="multi_riq","RaoQ",overlay=T)

visreg(ID_Plant[[5]],by="multi_riq","RaoQ",overlay=T,breaks=c(24,480))

visreg(ID_Plant[[3]],cond=list(multi_riq=100),"FRic.ab",overlay=T)#definir valor para o tamanho da rede

summary(ID_Plant[[1]])

visreg(ID_Plant[[2]],"FDis.ab")

visreg(ID_Plant[[2]],by="multi_riq","FDis.ab",overlay=T)

visreg(ID_Plant[[2]],by="multi_riq","FDis.ab",overlay=T,breaks=c(24,480))

visreg(ID_Plant[[2]],cond=list(multi_riq=100),"FDis.ab",overlay=T)#definir valor para o tamanho da rede


summary(data$multi_riq)

####SÓ ABUNDANCIA

ID_Plant<-list()
ID_Plant[[1]]<-lm(ID ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
ID_Plant[[2]]<-lm(ID ~ FDis.ab+multi_riq, data=data)
#ID_Plant[[2]]<-lm(ID ~ FRic.ab+multi_riq,data=data)
ID_Plant[[3]]<-lm(ID ~ FEve.ab+multi_riq, data=data)
#ID_Plant[[5]]<-lm(ID ~ RaoQ+multi_riq, data=data)
#ID_Plant[[4]]<-lm(ID ~ FDiv.ab+multi_riq, data=data)

ID_Plant[[4]]<-lm(ID ~ multi_riq, data=data)

head(model.sel(ID_Plant, rank="AICc"),2)
summary(ID_Plant[[4]])
model.sel(ID_Plant, rank="AICc")
visreg(ID_Plant[[2]],"FDis.ab")

visreg(ID_Plant[[2]],by="multi_riq","FDis.ab",overlay=T)

visreg(ID_Plant[[2]],by="multi_riq","FDis.ab",overlay=T,breaks=c(24,480))

p<-visreg(ID_Plant[[2]],cond=list(multi_riq=100),"FDis.ab",ylab="Diversidade de interações (ID)",xlab="Dispersão funcional (FDis) das plantas ponderada pela abundância",overlay=T)#definir valor para o tamanho da rede
library(ggplot2)  # Se ainda não carregado

# Extrair o gráfico do objeto visreg
meu_plot <- p$fit

# Plotar o gráfico
plot <- ggplot(data = meu_plot, aes(x = FDis.ab, y = ID)) +
  geom_line(aes(y = visregFit), color = "blue") +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = "blue", alpha = 0.2) +
  labs(x = "Dispersão funcional (FDis) das plantas ponderada pela abundância",
       y = "Diversidade de interações (ID)") +
  theme_minimal()+
  ggtitle("(a)")

plot

#SÓ PRESENÇA E AUSENCIA

ID_Plant<-list()
ID_Plant[[1]]<-lm(ID ~ 1, data=data)
#ID_Plant[[2]]<-lm(ID ~ FRic.pa+multi_riq, data=data)
ID_Plant[[2]]<-lm(ID ~ FDis.pa+multi_riq, data=data)
ID_Plant[[3]]<-lm(ID ~ FEve.pa+multi_riq, data=data)
#ID_Plant[[5]]<-lm(ID ~ FDiv.pa+multi_riq, data=data)
ID_Plant[[4]]<-lm(ID ~ multi_riq, data=data)

summary(ID_Plant[[4]])
model.sel(ID_Plant, rank="AICc")

visreg(ID_Plant[[2]],"FDis.pa")

visreg(ID_Plant[[2]],by="multi_riq","FDis.pa",overlay=T)

visreg(ID_Plant[[2]],by="multi_riq","FDis.pa",overlay=T,breaks=c(24,480))

p1<-visreg(ID_Plant[[2]],cond=list(multi_riq=100),"FDis.pa",overlay=T, ylab="Diversidade de interações (ID)",xlab="Dispersão funcional (FDis) das plantas por presença e ausência")
# Extrair o gráfico do objeto visreg
plot_p1 <- p1$fit

# Plotar o gráfico
plot1 <- ggplot(data = plot_p1, aes(x = FDis.pa, y = ID)) +
  geom_line(aes(y = visregFit), color = "blue") +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = "blue", alpha = 0.2) +
  labs(x = "Dispersão funcional (FDis) das plantas por presença e ausência",
       y = "Diversidade de interações (ID)") +
  theme_minimal()+
  ggtitle("(b)")
plot1

library(gridExtra)
dois<-grid.arrange(plot, plot1, ncol = 2)

ggsave("fdis_planta_pa_ab.png", plot = dois, width = 10, height = 6, dpi = 300)

###plot NSIZE
p3<-visreg(ID_Plant[[4]],cond=list(multi_riq=100),"multi_riq",overlay=T, ylab="Diversidade de interações (ID)",xlab="Tamanho da rede (NSize)")
# Extrair o gráfico do objeto visreg
plot_p1 <- p3$fit

# Plotar o gráfico
plot3 <- ggplot(data = plot_p1, aes(x = multi_riq, y = ID)) +
  geom_line(aes(y = visregFit), color = "blue") +
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), fill = "blue", alpha = 0.2) +
  labs(x = "Tamanho da rede (NSize)",
       y = "Diversidade de interações (ID)") +
  theme_minimal()

plot3


ggsave("nsize_ggplot.png", plot = plot3, width = 10, height = 6, dpi = 300)



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

#

IE_sum_Plant<-list()
IE_sum_Plant[[1]]<-lm(IE_sum ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
IE_sum_Plant[[2]]<-lm(IE_sum ~ FDis.ab+multi_riq, data=data)
IE_sum_Plant[[3]]<-lm(IE_sum ~ FRic.ab+multi_riq,data=data)
IE_sum_Plant[[4]]<-lm(IE_sum ~ FEve.ab+multi_riq, data=data)
IE_sum_Plant[[5]]<-lm(IE_sum ~ RaoQ+multi_riq, data=data)
IE_sum_Plant[[6]]<-lm(IE_sum ~ FDiv.ab+multi_riq, data=data)
IE_sum_Plant[[7]]<-lm(IE_sum ~ FRic.pa+multi_riq, data=data)
IE_sum_Plant[[8]]<-lm(IE_sum ~ FDis.pa+multi_riq, data=data)
IE_sum_Plant[[9]]<-lm(IE_sum ~ FEve.pa+multi_riq, data=data)
IE_sum_Plant[[10]]<-lm(IE_sum ~ FDiv.pa+multi_riq, data=data)
IE_sum_Plant[[11]]<-lm(IE_sum ~ multi_riq, data=data)

###IE ABUNDANCIA
IE_sum_Plant<-list()
IE_sum_Plant[[1]]<-lm(IE_sum ~ 1, data=data)
##Richness and abundance of data flora, those pollinated or seed dispersed by animals
IE_sum_Plant[[2]]<-lm(IE_sum ~ FDis.ab+multi_riq, data=data)

IE_sum_Plant[[3]]<-lm(IE_sum ~ FEve.ab+multi_riq, data=data)

IE_sum_Plant[[4]]<-lm(IE_sum ~ multi_riq, data=data)
model.sel(IE_sum_Plant, rank="AICc") #Best: NULL

summary(IE_sum_Plant[[4]])

###IE PRESENÇA E AUSENCIA
IE_sum_Plant<-list()
IE_sum_Plant[[1]]<-lm(IE_sum ~ 1, data=data)
IE_sum_Plant[[2]]<-lm(IE_sum ~ FDis.pa+multi_riq, data=data)
IE_sum_Plant[[3]]<-lm(IE_sum ~ FEve.pa+multi_riq, data=data)
IE_sum_Plant[[4]]<-lm(IE_sum ~ multi_riq, data=data)
model.sel(IE_sum_Plant, rank="AICc") #Best: NULL

summary(IE_sum_Plant[[2]])

car::Anova(IE_sum_Plant[[4]])

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(IE_sum_Plant, rank="AICc") #Best: NULL

summary(IE_sum_Plant[[1]])

visreg(IE_sum_Plant[[4]],"multi_riq")

visreg(IE_sum_Plant[[4]],by="multi_riq","multi_riq",overlay=T)

visreg(IE_sum_Plant[[4]],by="multi_riq","multi_riq",overlay=T,breaks=c(24,480))

visreg(IE_sum_Plant[[4]],cond=list(multi_riq=100),"multi_riq",overlay=T,xlab="Tamanho da rede (Nsize)",ylab="Homogeneidade de interações (IE)")#definir valor para o tamanho da rede
