library(vegan)

read.table("dados_moran.txt", h=T)->var_moran
xy<-var_moran[,1:2]

xy.dists <- as.matrix(dist(xy))
dim(xy.dists)
xy.dists.inv <- 1/xy.dists
diag(xy.dists.inv) <- 0

xy.dists.inv[1:14, 1:14]

library(ape)
names(var_moran)
Moran.I(var_moran$w_con, xy.dists.inv)#w_connectance
Moran.I(var_moran$ID, xy.dists.inv)#Shannon
Moran.I(var_moran$IE, xy.dists.inv)#IE
Moran.I(var_moran$H2, xy.dists.inv)#H2

#Sem autoco
