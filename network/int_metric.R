
library(bipartite)


MyData <- read.csv(file="Re_int_2024.csv", header=TRUE, sep=",")
MyData

web <- frame2webs(MyData, varnames = c("lower", "higher", "webID"), type.out = "list", emptylist = TRUE)
web


#######networklevel############
redeAll<-networklevel(web[["CN14"]], index="ALLBUTDD", level="both", weighted=TRUE,
                      ISAmethod="Bluethgen", SAmethod = "Bluethgen", extinctmethod = "r",
                      nrep = 1000, CCfun=median, dist="horn", normalise=TRUE, empty.web=TRUE,
                      logbase="e", intereven="sum", H2_integer=TRUE, fcweighted=TRUE,
                      fcdist="euclidean", legacy=FALSE)


write.csv(redeAll, "netlev_CN14.csv")


#######specieslevel#########

MyData <- read.csv(file="Re_int2_2024.csv", header=TRUE, sep=",")
MyData


web_sp <- frame2webs(MyData, varnames = c("lower", "higher", "webID"), type.out = "list", emptylist = TRUE)
web_sp

spredeAll <- specieslevel(web_sp[["CN"]], index="d", level="lower", logbase=exp(1), low.abun=NULL,
                          high.abun=NULL, PDI.normalise=TRUE, PSI.beta=c(1,0), nested.method="NODF",
                          nested.normalised=TRUE, nested.weighted=TRUE, empty.web=TRUE)


write.csv(spredeAll, "d_CN_low.csv")


#######plot#########

install.packages("bipartiteD3")
library(bipartiteD3)
library(bipartite)

MyData <- read.csv(file="Re_int2_2024.csv", header=TRUE, sep=",")
MyData

web_splot<- frame2webs(MyData, varnames = c("lower", "higher", "webID"), type.out = "list", emptylist = TRUE)
web_splot

bipartite_D3(web_splot)

plotweb(sortweb(web_splot[["CN"]], sort.order="dec", sequence=NULL), y.width.low = 0.3,
        y.width.high = 0.3, arrow="up", col.high = "black",
        col.low="black", bor.col.interaction ="black", method="normal", text.rot=90)





