#code for fig.2
library(ggplot2)
library(nlme)
library(colorRamps)
library(reshape2)
library(scales)
library(plotrix)
library(patchwork)
data<-read.table("Data1-plot.csv",header = TRUE, sep = ",")
data$year<-factor(data$year)
top_bar1 <- function(ANPP){  
  return(mean(ANPP)+std.error(ANPP)) 
}
bottom_bar1 <- function(ANPP){
  return(mean(ANPP)-std.error(ANPP))
}
top_bar2 <- function(BNPP){  
  return(mean(BNPP)+std.error(BNPP)) 
}
bottom_bar2 <- function(BNPP){
  return(mean(BNPP)-std.error(BNPP))
}
top_bar3 <- function(PAR){  
  return(mean(PAR)+std.error(PAR)) 
}
bottom_bar3 <- function(PAR){
  return(mean(PAR)-std.error(PAR))
}
plot1<-ggplot(data=data,aes(x=NR,y=ANPP,color=year,fill=year))+
  facet_wrap(~np,scales="free")+
  stat_summary(geom = 'point',size=3)+
  stat_summary(geom = 'errorbar',width=.05,cex=0.5,fun.min = bottom_bar1,fun.max = top_bar1)+
  geom_smooth(data=subset(subset(data,np=="N1P1"),year=="2019"),method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Above-ground biomass (g)'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:600*200)+coord_cartesian(ylim = c(0, 600))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot1
plot2<-ggplot(data=data,aes(x=NR,y=BNPP,color=year,fill=year))+
  facet_wrap(~np,scales="free")+
  stat_summary(geom = 'point',size=3)+
  stat_summary(geom = 'errorbar',width=.05,cex=0.5,fun.min = bottom_bar2,fun.max = top_bar2)+
  geom_smooth(data=subset(subset(data,np=="N1P1"),year=="2018"),method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Below-ground biomass (g)'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:400*100)+coord_cartesian(ylim = c(0, 400))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot2
plot3<-ggplot(data=data,aes(x=NR,y=PAR,color=year,fill=year))+
  facet_wrap(~np,scales="free")+
  stat_summary(geom = 'point',size=3)+
  stat_summary(geom = 'errorbar',width=.05,cex=0.5,fun.min = bottom_bar3,fun.max = top_bar3)+
  geom_smooth(data=subset(subset(data,np==c("N1P0","N1P1")),year=="2019"),method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('Fraction of PAR to the surface (%)'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:70*20)+coord_cartesian(ylim = c(0, 70))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot3
plot1+plot2+plot3+ plot_layout(ncol = 1, byrow = FALSE)

#code for fig.3
library(ggplot2)
library(nlme)
library(colorRamps)
library(reshape2)
library(scales)
library(plotrix)
library(patchwork)
data1<-read.table("Data1-plot.csv",header = TRUE, sep = ",")
data1$year<-factor(data1$year)
top_bar <- function(ESNpie){  
  return(mean(ESNpie)+std.error(ESNpie)) 
}
bottom_bar <- function(ESNpie){
  return(mean(ESNpie)-std.error(ESNpie))
}
plot1<-ggplot(data=data1,aes(x=NR,y=ESNpie,color=year,fill=year))+
  facet_wrap(~np,scales="free")+
  stat_summary(geom = 'point',size=3)+
  stat_summary(geom = 'errorbar',width=.05,cex=0.5,fun.min = bottom_bar,fun.max = top_bar)+
  geom_smooth(data=subset(subset(data,np=="N0P0"),year=="2020"),method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('¦Á diversity'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:4*1)+coord_cartesian(ylim = c(1, 4))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot1
data2<-read.table("Data2-treatment.csv",header = TRUE, sep = ",")
data2$year<-factor(data2$year)
plot2<-ggplot(data=data2,aes(x=NR,y=beta,color=year,fill=year))+
  facet_wrap(~np,scales="free")+geom_point(size=3)+
  geom_smooth(data=rbind(subset(subset(data,np=="N0P0"),year=="2018"),
                         subset(subset(data,year==c("2018","2019")),np=="N1P0"),
                         subset(data,np=="N1P1")),
              method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('¦Â diversity'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:10*0.1)+coord_cartesian(ylim = c(0.1, 0.5))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=16,face="bold"))
plot2
plot3<-ggplot(data=data2,aes(x=NR,y=gamma,color=year,fill=year))+
  facet_wrap(~np,scales="free")+geom_point(size=3)+
  geom_smooth(data=rbind(subset(subset(data,np=="N0P0"),year=="2018"),
                         subset(subset(data,year==c("2018","2020")),np=="N1P0"),
                         subset(subset(data,year==c("2018","2020")),np=="N1P1")),
              method = 'lm',se = TRUE, level=0.95)+
  scale_colour_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  scale_fill_manual( values = c("#0072b2", "#009e73","#e69f00"))+
  xlab(bquote('Number of added nutrients except N & P')) +ylab(bquote('¦Ã diversity'))+
  scale_x_continuous(limits = c(0,6), breaks=0:6*1)+
  scale_y_continuous(breaks=0:45*5)+coord_cartesian(ylim = c(20, 45))+
  theme_bw()+
  theme(text=element_text(size=20,  family="sans", face="bold"))+
  theme(legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size=16, face="bold"),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot3
plot1+plot2+plot3+ plot_layout(ncol = 1, byrow = FALSE)

#SEM code for Fig.4
library(piecewiseSEM)
library(nlme)
data<-read.table("Data1-plot.csv",header = TRUE, sep = ",")
data$block<-factor(data$block)
data$ANPP<-log(data$ANPP)
data$PAR<-(data$PAR)/100
model1<-psem(lme(ANPP~NP+BNPP+NR,random=~1|block,data=data),
             lme(BNPP~NP+NR,random=~1|block,data=data),
             lme(PAR~ANPP+NP+NR,random=~1|block,data=data),
             lme(ESNpie~NP+PAR+ANPP+NR,random=~1|block,data=data))
summary(model1)
model1.dsep<-dSep(model1,conditioning=T)
model2<-psem(lme(ANPP~NP+BNPP+NR,random=~1|block,data=data),
             lme(BNPP~NP+NR,random=~1|block,data=data),
             lme(PAR~ANPP+NP+NR,random=~1|block,data=data),
             lm(beta~NP+PAR+ANPP+NR,data=data))
summary(model2)
model2.dsep<-dSep(model2,conditioning=T)
model3<-psem(lme(ANPP~NP+BNPP+NR,random=~1|block,data=data),
             lme(BNPP~NP+NR,random=~1|block,data=data),
             lme(PAR~ANPP+NP+NR,random=~1|block,data=data),
             lm(gamma~NP+PAR+ANPP+NR,data=data))
summary(model3)
model3.dsep<-dSep(model3,conditioning=T)

#lm code for Tables 1 & 2
library(nlme)
library(MuMIn)
data1<-read.table("Data1-plot.csv",header = TRUE, sep = ",")
data1$year<-factor(data1$year)
data1$block<-factor(data1$block)
data1$NP<-factor(data1$NP)
data2<-read.table("Data2-treatment.csv",header = TRUE, sep = ",")
data2$year<-factor(data2$year)
data2$NP<-factor(data2$NP)
fullmodel1<-lme(ANPP~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year/block,data1)
summary(fullmodel1)
model_selection1 <- dredge(fullmodel1, evaluate = TRUE, rank = "AICc",REML=F)
model_selection1
subset(model_selection1, delta<2)
best_est1<-model.avg(model_selection1,subset= delta <2, revised.var = TRUE)
summary(best_est1)
importance(best_est1)
fullmodel2<-lme(BNPP~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year/block,data1)
summary(fullmodel2)
model_selection2 <- dredge(fullmodel2, evaluate = TRUE, rank = "AICc",REML=F)
model_selection2
subset(model_selection2, delta<2)
best_est2<-model.avg(model_selection2, revised.var = TRUE)
summary(best_est2)
importance(best_est2)
fullmodel3<-lme(PAR~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year/block,data1)
summary(fullmodel3)
model_selection3 <- dredge(fullmodel3, evaluate = TRUE, rank = "AICc",REML=F)
model_selection3
subset(model_selection3, delta<2)
best_est3<-model.avg(model_selection3,subset= delta <2, revised.var = TRUE)
summary(best_est3)
importance(best_est3)
fullmodel4<-lme(ESNpie~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year/block,data1)
summary(fullmodel4)
model_selection4 <- dredge(fullmodel4, evaluate = TRUE, rank = "AICc",REML=F)
model_selection4
subset(model_selection4, delta<2)
best_est4<-model.avg(model_selection4,subset= delta <2, revised.var = TRUE)
summary(best_est4)
importance(best_est4)
fullmodel5<-lme(beta~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year,data2)
summary(fullmodel5)
model_selection5 <- dredge(fullmodel5, evaluate = TRUE, rank = "AICc",REML=F)
model_selection5
subset(model_selection5, delta<2)
best_est5<-model.avg(model_selection5, revised.var = TRUE)
summary(best_est5)
importance(best_est5)
fullmodel6<-lme(gamma~NP+NR+NP:NR+NP:duration+NR:duration+NP:NR:duration, random=~1|year,data2)
summary(fullmodel6)
model_selection6 <- dredge(fullmodel6, evaluate = TRUE, rank = "AICc",REML=F)
model_selection6
subset(model_selection6, delta<2)
best_est6<-model.avg(model_selection6,subset= delta <2, revised.var = TRUE)
summary(best_est6)
importance(best_est6)
fit<-glm(gamma~NP+NR+NP:NR+NP:duration+NR:duration,family=quasipoisson,data=data2)
summary(fit)
exp(coef(fit))