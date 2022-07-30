#MCMCglmm code for Table1
library(ape)
library(MCMCglmm)
library(coda)
library(nlme)
data <- read.csv("Data1-species.csv", h=T)
names(data)  
head(data)
str(data)
cols <- c( 1:5)
data[cols] <- lapply(data[cols], as.factor)
str(data)
aggregate(data$Plots, by=list(data$Block, data$NAR), length)
nrow(data)
length(levels(as.factor(data$Plots)))
data$pres01 <- ceiling(data$presence)
str(data$pres01)
data$pres01  <- as.factor(data$pres01)
head(data)
tail(data)
data$logGS <- log(data$Cvalue)
tree <-read.nexus("tree_darwincn_proportional70.nex")
inv.phylo <- inverseA(tree,  scale=TRUE)
colnames(data)[1] <- "animal"
prior_idhBin <- list(R = list(V = 1, fix=1),
                     G = list(G1 = list(nu = 0.002, V = 1),
                              G2 = list(V = diag(5),
                                        nu = 1, alpha.mu=rep(0, 5),
                                        alpha.V=diag(5))))
prior_idhBin_b <- list(R = list(V = 1, fix=1),
                       G = list(G1 = list(nu = 0.002, V = 1),
                                G2 = list(V = diag(5),
                                          nu = 1, alpha.mu=rep(0, 5),
                                          alpha.V=diag(5)),
                                G3 = list(nu = 0.002, V = 1)))
set.seed(4222)
mbin3d <- MCMCglmm(pres01 ~ NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=data,
                   random= ~animal + idh(Block):species,
                   ginverse=list(animal= inv.phylo$Ainv),
                   verbose=T,  prior= prior_idhBin, family = "ordinal",
                   nitt=1000000, thin= 10, burnin=100000) 
summary(mbin3d)
set.seed(4321)
mbin3db <- MCMCglmm(pres01 ~ NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=data,
                    random= ~animal + idh(Block):species,
                    ginverse=list(animal= inv.phylo$Ainv),
                    verbose=T,  prior= prior_idhBin, family = "ordinal",
                    nitt=1000000, thin= 10, burnin=100000)
set.seed(1234)
mbin3dc <- MCMCglmm(pres01 ~ NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=data,
                    random= ~animal + idh(Block):species,
                    ginverse=list(animal= inv.phylo$Ainv),
                    verbose=T,  prior= prior_idhBin, family = "ordinal",
                    nitt=1000000, thin= 10, burnin=100000)
set.seed(2444)
mbin3dd <- MCMCglmm(pres01 ~ NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=data,
                    random= ~animal + idh(Block):species,
                    ginverse=list(animal= inv.phylo$Ainv),
                    verbose=T,  prior= prior_idhBin, family = "ordinal",
                    nitt=1000000, thin= 10, burnin=100000)
chainList_binModel <- mcmc.list(mbin3d $Sol, mbin3db $Sol, mbin3dc $Sol, mbin3dd$Sol)
gelman.diag(chainList_binModel)
plot(mcmc.list(mbin3db $Sol, mbin3dc $Sol))
plot(mcmc.list(mbin3dc $Sol, mbin3dd $Sol))
plot(mcmc.list(mbin3db $VCV, mbin3dc $VCV))
plot(mcmc.list(mbin3dc $VCV, mbin3dd $VCV))
#MCMCglmm code for Table2
totbm <- aggregate(data$Biomass, by=list(data$Plots), sum)
head(totbm)
colnames(totbm) <- c("Plots", "totBM")
data <- merge(data, totbm, by="Plots")
head(data)
data$bmpc <- data$Biomass/data$totBM *100
hist(data$bmpc)
sub10 <- subset(data, bmpc >= 10)
length(levels(factor(sub10$species)))
totspp <- aggregate(as.numeric(as.character(data$pres01)), by=list(data$Plots), sum)
head(totspp)
colnames(totspp) <- c("Plots", "totspp")
data <- merge(data, totspp, by="Plots")
keep_pc10 <- levels(factor(sub10 $species))
range(sub10$Cvalue)
dat <- data[data$species %in% keep_pc10,]
nrow(dat)
plot(dat $Biomass ~ jitter(dat $NAR), cex=0.7)
abline(lm(dat $Biomass ~(dat $NAR)))
plot(dat $Biomass ~ dat $logGS, cex=0.7)
plot(dat $logGS ~ jitter(dat $NAR), cex=0.7)
str(dat)
dat$pres01<-as.numeric(dat$pres01)
totspp2 <- aggregate(dat$pres01, by = list(dat$Plot), sum)
colnames(totspp2) <- c("Plots", "n_spp")
dat <- merge(dat, totspp2)
head(dat)
range(dat$n_spp)
boxplot(n_spp ~ NAR, data = dat)
plot(jitter(n_spp) ~ jitter(NAR), data = dat)
abline(lm(n_spp ~ NAR, data=dat))
prior_idh2 <- list(R = list(V = diag(5), nu = 0.002),
                   G = list(G1 = list(nu = 0.002, V = 1),
                            G2 = list(V = diag(5),
                                      nu = 0.002,
                                      alpha.mu=rep(0, 5),
                                      alpha.V=diag(5))))
set.seed(4222)
m2f <- MCMCglmm(Biomass ~  NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=dat,
                random= ~animal + idh(Block):species,
                ginverse=list(animal= inv.phylo$Ainv),
                verbose=T,  rcov= ~idh(Block):units,  prior= prior_idh2,
                nitt= 1000000, thin= 10, burnin= 100000)
summary(m2f)
plot(m2f $VCV)
plot(m2f $Sol)
autocorr.diag(m2f $Sol)
autocorr.diag(m2f $VCV)
heidel.diag(m2f $VCV)
heidel.diag(m2f $Sol)
geweke.plot(m2f $VCV)
geweke.plot(m2f $Sol)
geweke.diag(m2f $VCV) 
geweke.diag(m2f $Sol)
set.seed(4321)
m2fb<- MCMCglmm(Biomass ~  NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=dat,
                 random= ~animal + idh(Block):species,
                 ginverse=list(animal= inv.phylo$Ainv),
                 verbose=T,  rcov= ~idh(Block):units,  prior= prior_idh2,
                 nitt= 1000000, thin= 10, burnin= 100000)
set.seed(1234)
m2fc <- MCMCglmm(Biomass ~  NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=dat,
                 random= ~animal + idh(Block):species,
                 ginverse=list(animal= inv.phylo$Ainv),
                 verbose=T,  rcov= ~idh(Block):units,  prior= prior_idh2,
                 nitt= 1000000, thin= 10, burnin= 100000)
set.seed(2444)
m2fd <- MCMCglmm(Biomass ~  NAR * log(Cvalue) + NP + NP : log(Cvalue),  data=dat,
                 random= ~animal + idh(Block):species,
                 ginverse=list(animal= inv.phylo$Ainv),
                 verbose=T,  rcov= ~idh(Block):units,  prior= prior_idh2,
                 nitt= 1000000, thin= 10, burnin= 100000)
chainList_binModel <- mcmc.list(m2f$Sol, m2fb$Sol, m2fc$Sol, m2fd$Sol)
gelman.diag(chainList_binModel)
plot(mcmc.list(m2fb $Sol, m2fc $Sol))
plot(mcmc.list(m2fc $Sol, m2fd $Sol))
plot(mcmc.list(m2fb $VCV, m2fc $VCV))
plot(mcmc.list(m2fc $VCV, m2fd $VCV))
#Plotting code for fig.1
library(ggplot2)
library(nlme)
library(colorRamps)
library(reshape2)
library(scales)
library(plotrix)
library(patchwork)
data<-read.table("Data2-community.csv",header = TRUE, sep = ",")
top_bar1 <- function(ANPP){  
  return(mean(ANPP)+std.error(ANPP)) 
}
bottom_bar1 <- function(ANPP){
  return(mean(ANPP)-std.error(ANPP))
}
top_bar2 <- function(LA5){  
  return(mean(LA5)+std.error(LA5)) 
}
bottom_bar2 <- function(LA5){
  return(mean(LA5)-std.error(LA5))
}
top_bar3 <- function(SA5){  
  return(mean(SA5)+std.error(SA5)) 
}
bottom_bar3 <- function(SA5){
  return(mean(SA5)-std.error(SA5))
}
plot1<-ggplot(data=data,aes(x=NR,y=ANPP,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar1,fun.max = top_bar1)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('ANPP ('*g~m ^ -2*')'))+
  annotate("text",x=0,y=500,parse=TRUE,label='a',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:500*100)+coord_cartesian(ylim = c(0, 500))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.position = "none",
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot1
plot2<-ggplot(data=data,aes(x=NR,y=LA5,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar2,fun.max = top_bar2)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('ANPP ('*g~m ^ -2*')'))+
  annotate("text",x=1.3,y=500,parse=TRUE,label='b~(large~GS)',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:500*100)+coord_cartesian(ylim = c(0, 500))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.position = "none",
        axis.text.x=element_text(size=16,face="bold"),axis.text.y=element_blank(),
        axis.title.x=element_text(size=16,face="bold"),axis.title.y=element_blank())
plot2
plot3<-ggplot(data=data,aes(x=NR,y=SA5,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar3,fun.max = top_bar3)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('ANPP ('*g~m ^ -2*')'))+
  annotate("text",x=1.3,y=500,parse=TRUE,label='c~(small~GS)',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:500*100)+coord_cartesian(ylim = c(0, 500))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title =element_blank(),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text.x=element_text(size=16,face="bold"),axis.text.y=element_blank(),
        axis.title=element_blank())
plot3
plot1+plot2+plot3+plot_layout(ncol = 3, byrow = FALSE)
#Plotting code for fig.2
top_bar1 <- function(Richness){  
  return(mean(Richness)+std.error(Richness)) 
}
bottom_bar1 <- function(Richness){
  return(mean(Richness)-std.error(Richness))
}
top_bar2 <- function(LR5){  
  return(mean(LR5)+std.error(LR5)) 
}
bottom_bar2 <- function(LR5){
  return(mean(LR5)-std.error(LR5))
}
top_bar3 <- function(SR5){  
  return(mean(SR5)+std.error(SR5)) 
}
bottom_bar3 <- function(SR5){
  return(mean(SR5)-std.error(SR5))
}
plot1<-ggplot(data=data,aes(x=NR,y=Richness,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar1,fun.max = top_bar1)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Plant species richness'))+
  annotate("text",x=0,y=25,parse=TRUE,label='a',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:25*5)+coord_cartesian(ylim = c(0, 25))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.position = "none",
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot1
plot2<-ggplot(data=data,aes(x=NR,y=LR5,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar2,fun.max = top_bar2)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Plant species richness'))+
  annotate("text",x=1.3,y=25,parse=TRUE,label='b~(large~GS)',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:25*5)+coord_cartesian(ylim = c(0, 25))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.position = "none",
        axis.text.x=element_text(size=16,face="bold"),axis.text.y=element_blank(),
        axis.title.x=element_text(size=16,face="bold"),axis.title.y=element_blank())
plot2
plot3<-ggplot(data=data,aes(x=NR,y=SR5,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar3,fun.max = top_bar3)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Plant species richness'))+
  annotate("text",x=1.3,y=25,parse=TRUE,label='c~(small~GS)',size=8)+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:25*5)+coord_cartesian(ylim = c(0, 25))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title =element_blank(),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text.x=element_text(size=16,face="bold"),axis.text.y=element_blank(),
        axis.title=element_blank())
plot3
plot1+plot2+plot3+plot_layout(ncol = 3, byrow = FALSE)
#Plotting code for fig.3
top_bar <- function(PAR){  
  return(mean(PAR)+std.error(PAR)) 
}
bottom_bar <- function(PAR){
  return(mean(PAR)-std.error(PAR))
}
plot<-ggplot(data=data,aes(x=NR,y=PAR,fill=np))+
  stat_summary(geom = 'errorbar',width=0.15,cex=0.6,fun.min = bottom_bar,fun.max = top_bar)+
  stat_summary(geom = 'point',size=5,color="black",shape=21)+
  scale_fill_manual( values = c("#c0c0c0", "#94ec7e","#067d55"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Fraction of PAR to the surface (%)'))+
  scale_x_continuous(breaks=0:6*1)+
  scale_y_continuous(breaks=0:50*10)+coord_cartesian(ylim = c(0, 50))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title =element_blank(),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot
#SEM code for Figure 4
library(lavaan)
model='
LA5~NP
SA5~LA5+NP
PAR~LA5+SA5+NP
LR5~LA5+NP+PAR
SR5~SA5+NP+PAR'
gs<-sem(model=model, data=data)
summary(gs,rsq=T,standardize=T)
summary(gs,modindices=T)
fitMeasures(gs)