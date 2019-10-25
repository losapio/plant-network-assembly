#### Script by Gianalberto Losapio
#### cite as: The assembly of a plant network in alpine vegetation
#### Journal of Vegetation Science
#### Losapio et al. 2018

library(nlme)
library(spatstat)
library(bipartite)
library(igraph)
library(MASS)
library(car)
library(vegan)
library(FactoMineR)

load("/Users/losapiog/g/pre18/phd/r/spatial_ppnet/jvs_20180717.RData")
# global networks, i.e. overall scales = 1 value across 75 cm

396
496
517
534
u.net  <- matrix(NA,19,19, dimnames=c(list(sp.cod),list(sp.cod)))
pu.net <- matrix(NA,19,19, dimnames=c(list(sp.cod),list(sp.cod)))

counter <-0

for(k in 1:nsp){
for(z in 1:nsp){
	if(k!=z){
	counter <- counter+1; print(counter)
		u.net[k,z]  <- test.indep.env[[19*(k-1)+z]][[2]]$statistic$u
		pu.net[k,z] <- test.indep.env[[19*(k-1)+z]][[2]]$p.value
}}}

# quantitative global networks: u = association strength
glonet.qnt <- u.net
glonet.qnt <- ifelse(pu.net <= 0.05, glonet.qnt,0)
glonet.qnt <- ifelse(is.na(glonet.qnt), 0,glonet.qnt)

# qualitative global networks: u = 1
glonet.qlv <- u.net
glonet.qlv <- ifelse(pu.net <= 0.05, 1,0)
glonet.qlv<- ifelse(is.na(glonet.qlv), 0,glonet.qlv)
g.glonet.qlv <- graph.adjacency(glonet.qlv)

plot(g.glonet.qlv)

##

deg.glonet <- degree(g.glonet.qlv)

trans.obs <- transitivity(g.glonet.qlv, "global", isolates="zero")

### null models

null.glonet.er <- rep(list(NA), 999)
for(i in 1:999){null.glonet.er[[i]] <- sample_gnm(19, 36, directed=T)}

null.df <- data.frame(simu=1:999, modul=0, trans=0)

for(i in 1:999){
	print(i)
	null.df$trans[i] <- transitivity(null.glonet.er[[i]], "global", isolates="zero")
}

null.deg.glonet <- rep(list(NA), 999)
for(i in 1:999){null.deg.glonet[[i]] <- degree(null.glonet.er[[i]])}

## media e sd
# deg
null.deg.media <- rep(NA, nsp)
null.deg.sd    <- rep(NA, nsp)

vettorenull <- rep(NA, 999)
for(j in 1:19){
	for(i in 1:999){
		vettorenull[i]<-null.deg.glonet[[i]][j]
		}
	null.deg.media[j]	<-mean(vettorenull)
	null.deg.sd[j]		<-sd(vettorenull)
}

z.deg <- rep(NA, nsp)
z.deg <- (deg.glonet-null.deg.media)/null.deg.sd

sptr$zdeg<-z.deg
z.deg <- z.deg[order(z.deg)]

# calculate p value

i.trans <- rep(0,999)

for(i in 1:999){
	if(trans.obs > null.df$trans[i]) i.trans[i] <- 1
}

p.trans <- 1-(sum(i.trans)/1000)

# degree
p.deg.h <- rep(0,19)
p.deg.l <- rep(0,19)

for(j in 1:19){
	i.deg <- rep(0,999)
for(i in 1:999){
	if(deg.glonet[j] > null.deg.glonet[[i]][j]) i.deg[i] <- 1
}

p.deg.h[j] <- 1-(sum(i.deg)/1000)

	i.deg <- rep(0,999)
for(i in 1:999){
	if(deg.glonet[j] < null.deg.glonet[[i]][j]) i.deg[i] <- 1
}

p.deg.l[j] <- 1-(sum(i.deg)/1000)

}

sptr$p.deg.h<-p.deg.h
sptr$p.deg.l<-p.deg.l

#
####
sptr.scale <- sptr
rownames(sptr.scale)<-sptr$sp
#sptr.scale[13,c("lma","leaves")]<-NA
sptr.scale <- sptr.scale[-13,]

mod.pca<-PCA(sptr.scale[,c("abund","cover","diam","height","lma","leaves", "zdeg")], quanti.sup=7, scale=T)

summary(mod.pca)
mod.pca$eig

dimdesc(mod.pca, 1:2)
mod.pca$var

sptr.scale$pc1<-mod.pca$ind$coord[,1]
sptr.scale$pc2<-mod.pca$ind$coord[,2]

mod.zdeg<- lm(zdeg ~ pc1 + pc2, data=sptr.scale)
anova(mod.zdeg)
summary(mod.zdeg)