library(readr)
data <- read_delim("func_models.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

head(data)

hist(data$ID)

summary(data$ID)

max(data$multi_riq)

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

fit.normal=fitdist(data$ID[-28],"norm")
fit.lognormal=fitdist(data$ID[-28],"lnorm")
fit.gamma=fitdist(data$ID[-28],"gamma")
par(mfrow=c(1,2))

#plotando as tres distribuições
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição

###mejor ajuste: normal
par(mfrow=c(1,1))



#ID ~ FRic (Bees + Plants)
names(data)

###USANDO modelos nulos (randomização da abundancia na FD)

ID_fric<-list()
ID_fric[[1]]<-lm(ID ~ 1, data=data)


ID_fric[[2]]<-lm(ID ~ FRic.n.bee+FRic.n.plant+multi_riq, data=data)
ID_fric[[3]]<-lm(ID ~ FRic.n.bee+multi_riq, data=data)
ID_fric[[4]]<-lm(ID ~ FRic.n.plant+multi_riq, data=data)
ID_fric[[5]]<-lm(ID ~ FRic.n.bee, data=data)
ID_fric[[6]]<-lm(ID ~ FRic.n.plant, data=data)
ID_fric[[7]]<-lm(ID ~ FRic.n.bee+FRic.n.plant, data=data)

ID_fric[[8]]<-lm(ID ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_fric, rank="AICc") #Best: 
summary(ID_fric[[3]])




#ID ~ FEve (Bees + Plants)
names(data)

###

ID_fEve<-list()
ID_fEve[[1]]<-lm(ID ~ 1, data=data)


ID_fEve[[2]]<-lm(ID ~ FEve.n.bee+FEve.n.plant+multi_riq, data=data)
ID_fEve[[3]]<-lm(ID ~ FEve.n.bee+multi_riq, data=data)
ID_fEve[[4]]<-lm(ID ~ FEve.n.plant+multi_riq, data=data)
ID_fEve[[5]]<-lm(ID ~ FEve.n.bee, data=data)
ID_fEve[[6]]<-lm(ID ~ FEve.n.plant, data=data)
ID_fEve[[7]]<-lm(ID ~ FEve.n.bee+FEve.n.plant, data=data)

ID_fEve[[8]]<-lm(ID ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_fEve, rank="AICc") #Best:
summary(ID_fEve[[2]])
vis_f2a<-visreg(ID_fEve[[2]],"FEve.n.plant")
vis_f2b<-visreg(ID_fEve[[2]],"multi_riq")
vis_f2a_bee<-visreg(ID_fEve[[2]],"FEve.n.bee")
summary(ID_fEve[[3]])
Anova(ID_fEve[[2]])


library(ggplot2)
library(ggeffects)
library(ggthemes)

dados_vis <- vis_f2a$fit
dados_residuais <- vis_f2a$res

# Criar o gráfico com ggplot2
Fig2a_plant<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis, aes(x = FEve.n.plant, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis, aes(x = FEve.n.plant, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais, aes(x = FEve.n.plant, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "",
    y = "Interaction Diversity (ID)"
  ) +
  theme(axis.title = element_text()
  )


dados_vis <- vis_f2a_bee$fit
dados_residuais <- vis_f2a_bee$res


Fig2a_bee<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis, aes(x = FEve.n.bee, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis, aes(x = FEve.n.bee, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais, aes(x = FEve.n.bee, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "",
    y = "Interaction Diversity (ID)"
  ) +
  theme(axis.title = element_text()
  )


dados_vis <- vis_f2b$fit
dados_residuais <- vis_f2b$res

Fig2b<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis, aes(x = multi_riq, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis, aes(x = multi_riq, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais, aes(x = multi_riq, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "Network size",
    y = "Interaction Diversity (ID)"
  ) +
  theme(axis.title = element_text()
  )




#ID ~ FDis (Bees + Plants)
names(data)

###

ID_FDis<-list()
ID_FDis[[1]]<-lm(ID ~ 1, data=data)


ID_FDis[[2]]<-lm(ID ~ FDis.n.bee+FDis.n.plant+multi_riq, data=data)
ID_FDis[[3]]<-lm(ID ~ FDis.n.bee+multi_riq, data=data)
ID_FDis[[4]]<-lm(ID ~ FDis.n.plant+multi_riq, data=data)
ID_FDis[[5]]<-lm(ID ~ FDis.n.bee, data=data)
ID_FDis[[6]]<-lm(ID ~ FDis.n.plant, data=data)
ID_FDis[[7]]<-lm(ID ~ FDis.n.bee+FDis.n.plant, data=data)

ID_FDis[[8]]<-lm(ID ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(ID_FDis, rank="AICc") #Best: 
summary(ID_FDis[[8]])




####para H2

#####Model selection?
# Distr de la var de respuesta
hist(data$H2)
data_clean<-na.omit(data)

require(fitdistrplus)
fit.normal=fitdist(data_clean$H2[-6],"norm")
fit.lognormal=fitdist(data_clean$H2[-6],"lnorm")
fit.gamma=fitdist(data_clean$H2[-6],"gamma")
par(mfrow=c(1,2))

#plotando as tres distribuições
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição
#SE EU QUISER USAR O LOG OU O GAMA?
summary(data_clean$H2)
###mejor ajuste: normal
par(mfrow=c(1,1))

#ID ~ FRic (Bees + Plants)
names(data)

###USANDO modelos nulos (randomização da abundancia na FD)

H2_fric<-list()
H2_fric[[1]]<-lm(H2 ~ 1, data=data)


H2_fric[[2]]<-lm(H2 ~ FRic.n.bee+FRic.n.plant+multi_riq, data=data)
H2_fric[[3]]<-lm(H2 ~ FRic.n.bee+multi_riq, data=data)
H2_fric[[4]]<-lm(H2 ~ FRic.n.plant+multi_riq, data=data)
H2_fric[[5]]<-lm(H2 ~ FRic.n.bee, data=data)
H2_fric[[6]]<-lm(H2 ~ FRic.n.plant, data=data)
H2_fric[[7]]<-lm(H2 ~ FRic.n.bee+FRic.n.plant, data=data)

H2_fric[[8]]<-lm(H2 ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(H2_fric, rank="AICc") #Best: 
summary(H2_fric[[3]])




#ID ~ FEve (Bees + Plants)
names(data)

###

H2_fEve<-list()
H2_fEve[[1]]<-lm(H2 ~ 1, data=data)


H2_fEve[[2]]<-lm(H2 ~ FEve.n.bee+FEve.n.plant+multi_riq, data=data)
H2_fEve[[3]]<-lm(H2 ~ FEve.n.bee+multi_riq, data=data)
H2_fEve[[4]]<-lm(H2 ~ FEve.n.plant+multi_riq, data=data)
H2_fEve[[5]]<-lm(H2 ~ FEve.n.bee, data=data)
H2_fEve[[6]]<-lm(H2 ~ FEve.n.plant, data=data)
H2_fEve[[7]]<-lm(H2 ~ FEve.n.bee+FEve.n.plant, data=data)

H2_fEve[[8]]<-lm(H2 ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(H2_fEve, rank="AICc") #Best: 
summary(H2_fEve[[5]])
vis_f2c<-visreg(H2_fEve[[3]],"FEve.n.bee")
Anova(H2_fEve[[5]])

dados_vis2c <- vis_f2c$fit
dados_residuais2c <- vis_f2c$res

Fig2c<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis2c, aes(x = FEve.n.bee, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis2c, aes(x = FEve.n.bee, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais2c, aes(x = FEve.n.bee, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "Functional Evenness - Bees",
    y = "Network Specialization (H2’)"
  ) +
  theme(axis.title = element_text()
  )






#H2 ~ FDis (Bees + Plants)
names(data)

###

H2_FDis<-list()
H2_FDis[[1]]<-lm(H2 ~ 1, data=data)


H2_FDis[[2]]<-lm(H2 ~ FDis.n.bee+FDis.n.plant+multi_riq, data=data)
H2_FDis[[3]]<-lm(H2 ~ FDis.n.bee+multi_riq, data=data)
H2_FDis[[4]]<-lm(H2 ~ FDis.n.plant+multi_riq, data=data)
H2_FDis[[5]]<-lm(H2 ~ FDis.n.bee, data=data)
H2_FDis[[6]]<-lm(H2 ~ FDis.n.plant, data=data)
H2_FDis[[7]]<-lm(H2 ~ FDis.n.bee+FDis.n.plant, data=data)

H2_FDis[[8]]<-lm(H2 ~ multi_riq, data=data)

###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(H2_FDis, rank="AICc") #Best: 
summary(H2_FDis[[6]])




########CONECTANCIA PONDERADA


head(data)

hist(data$w_conec)

summary(data$w_conec)

require(fitdistrplus)
fit.normal=fitdist(data$w_conec,"norm")
fit.lognormal=fitdist(data$w_conec,"lnorm")
fit.gamma=fitdist(data$w_conec,"gamma")
par(mfrow=c(1,2))
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição


#W_conec ~ FRic
w_conec.Fric<-list()
w_conec.Fric[[1]]<-lm(w_conec ~ 1, data=data)


w_conec.Fric[[2]]<-lm(w_conec ~ FRic.n.bee+FRic.n.plant+multi_riq, data=data)
w_conec.Fric[[3]]<-lm(w_conec ~ FRic.n.bee+multi_riq, data=data)
w_conec.Fric[[4]]<-lm(w_conec ~ FRic.n.plant+multi_riq, data=data)
w_conec.Fric[[5]]<-lm(w_conec ~ FRic.n.bee+FRic.n.plant, data=data)
w_conec.Fric[[6]]<-lm(w_conec ~ FRic.n.bee, data=data)
w_conec.Fric[[7]]<-lm(w_conec ~ FRic.n.plant, data=data)

w_conec.Fric[[8]]<-lm(w_conec ~ multi_riq, data=data)



###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(w_conec.Fric, rank="AICc") #Best:

summary(w_conec.Fric[[8]])
visreg(w_conec.Fric[[8]],"multi_riq")



# W_conec ~ FEve
w_conec.FEve<-list()
w_conec.FEve[[1]]<-lm(w_conec ~ 1, data=data)


w_conec.FEve[[2]]<-lm(w_conec ~ FEve.n.bee+FEve.n.plant+multi_riq, data=data)
w_conec.FEve[[3]]<-lm(w_conec ~ FEve.n.bee+multi_riq, data=data)
w_conec.FEve[[4]]<-lm(w_conec ~ FEve.n.plant+multi_riq, data=data)
w_conec.FEve[[5]]<-lm(w_conec ~ FEve.n.bee+FEve.n.plant, data=data)
w_conec.FEve[[6]]<-lm(w_conec ~ FEve.n.bee, data=data)
w_conec.FEve[[7]]<-lm(w_conec ~ FEve.n.plant, data=data)

w_conec.FEve[[8]]<-lm(w_conec ~ multi_riq, data=data)



###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(w_conec.FEve, rank="AICc") #Best: 

summary(w_conec.FEve[[4]])
vis_f2e<-visreg(w_conec.FEve[[4]],"FEve.n.plant")
Anova(w_conec.FEve[[4]])

vis_f2f<-visreg(w_conec.FEve[[4]],"multi_riq")

dados_vis2e <- vis_f2e$fit
dados_residuais2e <- vis_f2e$res

Fig2e<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis2e, aes(x = FEve.n.plant, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis2e, aes(x = FEve.n.plant, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais2e, aes(x = FEve.n.plant, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "Functional Evenness - Plants",
    y = "Weighted Connectance (C)"
  ) +
  theme(axis.title = element_text()
  )


dados_vis2f <- vis_f2f$fit
dados_residuais2f <- vis_f2f$res

Fig2f<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis2f, aes(x = multi_riq, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis2f, aes(x = multi_riq, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais2f, aes(x = multi_riq, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "Network size",
    y = "Weighted Connectance (C)"
  ) +
  theme(axis.title = element_text()
  )


# W_conec ~ FDis
w_conec.FDis<-list()
w_conec.FDis[[1]]<-lmer(w_conec ~ 1, data=data)


w_conec.FDis[[2]]<-lm(w_conec ~ FDis.n.bee+FDis.n.plant+multi_riq, data=data)
w_conec.FDis[[3]]<-lm(w_conec ~ FDis.n.bee+multi_riq, data=data)
w_conec.FDis[[4]]<-lm(w_conec ~ FDis.n.plant+multi_riq, data=data)
w_conec.FDis[[5]]<-lm(w_conec ~ FDis.n.bee+FDis.n.plant, data=data)
w_conec.FDis[[6]]<-lm(w_conec ~ FDis.n.bee, data=data)
w_conec.FDis[[7]]<-lm(w_conec ~ FDis.n.plant, data=data)

w_conec.FDis[[8]]<-lm(w_conec ~ multi_riq, data=data)



###DBH for data plants, and for those pollinated or dispersed by animals:

model.sel(w_conec.FDis, rank="AICc") #Best: 

summary(w_conec.FDis[[8]])
visreg(w_conec.FDis[[8]],"multi_riq")
Anova(w_conec.FDis[[8]])



########### REGULARIDADE INTERAÇÔES

hist(data$IE_sum)

summary(data$IE_sum)


require(fitdistrplus)
fit.normal=fitdist(data$IE_sum[-28],"norm")
fit.lognormal=fitdist(data$IE_sum[-28],"lnorm")
fit.gamma=fitdist(data$IE_sum[-28],"gamma")
par(mfrow=c(1,2))
cdfcomp(list(fit.normal,fit.lognormal, fit.gamma),horizontals=F, addlegend=T,legendtext=c("Normal","log","gamma"))
qqcomp(list(fit.normal,fit.lognormal, fit.gamma),addlegend=T,legendtext=c("Normal","log", "gamma"))
gofstat(list(fit.normal, fit.lognormal, fit.gamma))$aic #o menor valor representa a melhor distribuição



#IE ~Fric
IE_FRic<-list()
IE_FRic[[1]]<-lm(IE ~ 1, data=data)


IE_FRic[[2]]<-lm(IE ~ FRic.n.bee+FRic.n.plant+multi_riq, data=data)
IE_FRic[[3]]<-lm(IE ~ FRic.n.bee+multi_riq, data=data)
IE_FRic[[4]]<-lm(IE ~ FRic.n.plant+multi_riq, data=data)
IE_FRic[[5]]<-lm(IE ~ FRic.n.bee+FRic.n.plant, data=data)
IE_FRic[[6]]<-lm(IE ~ FRic.n.bee, data=data)
IE_FRic[[7]]<-lm(IE ~ FRic.n.plant, data=data)

IE_FRic[[8]]<-lm(IE ~ multi_riq, data=data)



###DBH 

model.sel(IE_FRic, rank="AICc") #Best: Nsize

summary(IE_FRic[[8]])




#IE ~FEve
IE_FEve<-list()
IE_FEve[[1]]<-lm(IE ~ 1, data=data)


IE_FEve[[2]]<-lm(IE ~ FEve.n.bee+FEve.n.plant+multi_riq, data=data)
IE_FEve[[3]]<-lm(IE ~ FEve.n.bee+multi_riq, data=data)
IE_FEve[[4]]<-lm(IE ~ FEve.n.plant+multi_riq, data=data)
IE_FEve[[5]]<-lm(IE ~ FEve.n.bee+FEve.n.plant, data=data)
IE_FEve[[6]]<-lm(IE ~ FEve.n.bee, data=data)
IE_FEve[[7]]<-lm(IE ~ FEve.n.plant, data=data)

IE_FEve[[8]]<-lm(IE ~ multi_riq, data=data)



###DBH 

model.sel(IE_FEve, rank="AICc") #Best: Nsize

summary(IE_FEve[[8]])


#IE ~FDis
IE_FDis<-list()
IE_FDis[[1]]<-lm(IE ~ 1, data=data)


IE_FDis[[2]]<-lm(IE ~ FDis.n.bee+FDis.n.plant+multi_riq, data=data)
IE_FDis[[3]]<-lm(IE ~ FDis.n.bee+multi_riq, data=data)
IE_FDis[[4]]<-lm(IE ~ FDis.n.plant+multi_riq, data=data)
IE_FDis[[5]]<-lm(IE ~ FDis.n.bee+FDis.n.plant, data=data)
IE_FDis[[6]]<-lm(IE ~ FDis.n.bee, data=data)
IE_FDis[[7]]<-lm(IE ~ FDis.n.plant, data=data)

IE_FDis[[8]]<-lm(IE ~ multi_riq, data=data)



###DBH 

model.sel(IE_FDis, rank="AICc") #Best: Nsize

summary(IE_FDis[[8]])
vis_f2d<-visreg(IE_FDis[[8]], "multi_riq")
Anova(IE_FDis[[8]])
dados_vis2d <- vis_f2d$fit
dados_residuais2d <- vis_f2d$res

Fig2d<-ggplot() +
  # Adicionar a linha ajustada
  geom_line(data = dados_vis2d, aes(x = multi_riq, y = visregFit), 
            color = "blue", size = 1) +
  # Adicionar o intervalo de confiança
  geom_ribbon(data = dados_vis2d, aes(x = multi_riq, ymin = visregLwr, ymax = visregUpr), 
              fill = "blue", alpha = 0.2) +
  # Adicionar os pontos residuais
  geom_point(data = dados_residuais2d, aes(x = multi_riq, y = visregRes), 
             color = "gray50", alpha = 0.7) +
  # Personalizar o tema e os rótulos
  theme_classic() +
  labs(
    title = "",
    x = "Network size",
    y = "Interaction Evenness (IE)"
  ) +
  theme(axis.title = element_text()
  )

library(cowplot)
Figura.2 <- plot_grid(
  Fig2a_bee,Fig2a_plant,Fig2c,Fig2e, # Gráficos a serem combinados
  labels = c("a", "b", "c", "d"), # Rótulos para os gráficos
  label_size = 14, # Tamanho da fonte dos rótulos
  ncol = 2              # Número de colunas (os gráficos estarão lado a lado)
)

ggsave("Figura_2.tiff", Figura.2, units="mm",
       height = 130, width= 180, dpi=600,
       compression = "lzw")

Figura.2.1 <- plot_grid(
  Fig2a,Fig2b,Fig2c,Fig2d,Fig2e, # Gráficos a serem combinados
  ncol = 2              # Número de colunas (os gráficos estarão lado a lado)
)

Final<-ggdraw(Figura.2.1)+
  draw_label("a", x = .46, y = .96, size = 14, fontface="bold") +
  draw_label("b", x = .96, y = .96, size = 14, fontface="bold") +
  draw_label("c", x = .46, y = .63, size = 14, fontface="bold") +
  draw_label("d", x = .96, y = .63, size = 14, fontface="bold") +
  draw_label("e", x = .46, y = .29, size = 14, fontface="bold") +
  draw_label("f", x = .96, y = .29, size = 14, fontface="bold")
  
