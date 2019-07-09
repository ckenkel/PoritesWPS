#devtools::install_github("baptiste/egg")


library(lme4)
library(ggplot2)
library(car)
library(multcomp)
library(gtools)
library(plyr)
library(quantreg)
library(calibrate)
library(MASS)
library(AICcmodavg)
library(e1071)
library(nlme)
library(lmmfit)
library(MCMCglmm)
library(labdsv)
library(vegan)
library(plotrix)
library(pgirmess)
library(gridExtra)
library(egg)

#############################

###############################

setwd("") #location of input files - fill in here

#######################
#Read in data
########################

dat=read.csv("Patchy_Porites_Feb2019_for_R.csv",stringsAsFactors=TRUE)

summary(dat)

dat$condition=factor(dat$condition,levels=c("healthy","WPS"))

###############plots

pd <- position_dodge(.2)
p5<-ggplot(dat,aes(factor(condition),daycal_CaCO3.cm2.min,fill=condition))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=factor(sid)),position=pd)+
	geom_line(aes(group=factor(sid),colour=factor(sid)),position=pd)+
	labs(colour="genet")+
	#scale_colour_manual(values=c(point="gray60"))+
	labs(subtitle="Calcification",y=expression(CaCO["3"]/cm["2"]/min))+
	xlab(NULL)+
	#theme(axis.text.x=element_text(size=20))+
	theme_classic()+
	scale_fill_manual(values=c("snow4","snow2"),guide=FALSE)+
	theme(axis.text=element_text(size=11))+
	theme(plot.subtitle=element_text(hjust=0.5))

p1<-ggplot(dat,aes(factor(condition),zooxdcells.cm2,fill=condition))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=factor(sid)),position=pd,show.legend=FALSE)+
	geom_line(aes(group=factor(sid),colour=factor(sid)),position=pd,show.legend=FALSE)+
	#scale_fill_manual(guide=FALSE)+
	labs(subtitle="Symbiont Density", y=expression(cells/cm["2"]))+
	xlab(NULL)+
	#theme(axis.text.x=element_text(size=20))+
	theme_classic()+
	scale_fill_manual(values=c("snow4","snow2"),guide=FALSE)+
	theme(axis.text=element_text(size=11))+
	theme(plot.subtitle=element_text(hjust=0.5))
	
p2<-ggplot(dat,aes(factor(condition),chla_ug.cm2,fill=condition))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=factor(sid)),position=pd,show.legend=FALSE)+
	geom_line(aes(group=factor(sid),colour=factor(sid)),position=pd,show.legend=FALSE)+
	labs(colour="genet")+
	#scale_colour_manual(values=c(point="gray60"))+
	labs(subtitle="Chlorophyll a", y=expression(ug/cm["2"]))+
	xlab(NULL)+
	#theme(axis.text.x=element_text(size=20))+
	theme_classic()+
	scale_fill_manual(values=c("snow4","snow2"),guide=FALSE)+
	theme(axis.text=element_text(size=11))+
	theme(plot.subtitle=element_text(hjust=0.5))

p3<-ggplot(dat,aes(factor(condition),GrossPS_mgO2.cm2.min,fill=condition))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=factor(sid)),position=pd,show.legend=FALSE)+
	geom_line(aes(group=factor(sid),colour=factor(sid)),position=pd,show.legend=FALSE)+
	labs(colour="genet")+
	#scale_colour_manual(values=c(point="gray60"))+
	labs(subtitle="Gross Photosynthesis", y=expression(mgO["2"]/cm["2"]/min))+
	xlab(NULL)+
	#theme(axis.text.x=element_text(size=20))+
	theme_classic()+
	scale_fill_manual(values=c("snow4","snow2"),guide=FALSE)+
	theme(axis.text=element_text(size=11))+
	theme(plot.subtitle=element_text(hjust=0.5))

p4<-ggplot(dat,aes(factor(condition),protein_mg.cm2,fill=condition))+
	geom_boxplot(outlier.shape=16, notch=FALSE)+
	geom_point(aes(colour=factor(sid)),position=pd,show.legend=FALSE)+
	geom_line(aes(group=factor(sid),colour=factor(sid)),position=pd,show.legend=FALSE)+
	labs(colour="genet")+
	#scale_colour_manual(values=c(point="gray60"))+
	labs(subtitle="Soluble Protein",y=expression(mg/cm["2"]))+
	xlab(NULL)+
	#theme(axis.text.x=element_text(size=20))+
	theme_classic()+
	scale_fill_manual(values=c("snow4","snow2"),guide=FALSE)+
	theme(axis.text=element_text(size=11))+
	theme(plot.subtitle=element_text(hjust=0.5))

ggarrange(p1,p2,p3,p4,p5,ncol=2,labels=c('a','b','c','d','e'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))

########statistics on physiology variables

model1 <- lme(zooxdcells.cm2~condition,data=dat,random= ~1|sid,method="REML",na.action=na.omit) #change to different traits

par(mfrow=c(1,2))
#Diagnostics to determine if models are any good to begin with
qqPlot(residuals(model1),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal

qqnorm(model1, ~ranef(., level=1))


#check homoscedascitity by plotting residuals against fitted values
plot(model1,which=1) #want random scatter; no apparent trendline


summary(model1)
anova(model1)

Fixed effects: zooxdcells.cm2 ~ condition 
                   Value Std.Error DF   t-value p-value
(Intercept)      1803310  255403.1  7  7.060645   2e-04
conditionbleach -1544798  285255.7  7 -5.415486   1e-03

#########################################################
model2 <- lme(chla_ug.cm2~condition,data=dat,random= ~1|sid,method="REML",na.action=na.omit) #change to different traits


#Diagnostics to determine if models are any good to begin with
qqPlot(residuals(model2),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal

qqnorm(model2, ~ranef(., level=1))


#check homoscedascitity by plotting residuals against fitted values
plot(model2,which=1) #want random scatter; no apparent trendline


summary(model2)
anova(model2)

Fixed effects: chla_ug.cm2 ~ condition 
                    Value Std.Error DF   t-value p-value
(Intercept)      21.58604  1.875528  7 11.509312   0e+00
conditionbleach -18.42922  2.652398  7 -6.948136   2e-04

#########################################################
model3 <- lme(GrossPS_mgO2.cm2.min~condition,data=dat,random= ~1|sid,method="REML",na.action=na.omit) #change to different traits


#Diagnostics to determine if models are any good to begin with
qqPlot(residuals(model3),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal

qqnorm(model3, ~ranef(., level=1))


#check homoscedascitity by plotting residuals against fitted values
plot(model3,which=1) #want random scatter; no apparent trendline


summary(model3)
anova(model3)

Fixed effects: GrossPS_mgO2.cm2.min ~ condition 
                     Value Std.Error DF   t-value p-value
(Intercept)      1.2017362 0.1252366  7  9.595729   0e+00
conditionbleach -0.8978292 0.1508831  7 -5.950494   6e-04

#########################################################
model4 <- lme(protein_mg.cm2~condition,data=dat,random= ~1|sid,method="REML",na.action=na.omit) #change to different traits


#Diagnostics to determine if models are any good to begin with
qqPlot(residuals(model4),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal

qqnorm(model4, ~ranef(., level=1))


#check homoscedascitity by plotting residuals against fitted values
plot(model4,which=1) #want random scatter; no apparent trendline


summary(model4)
anova(model4)

Fixed effects: protein_mg.cm2 ~ condition 
                   Value Std.Error DF   t-value p-value
(Intercept)      623.639  75.08604  7  8.305658  0.0001
conditionbleach -256.745 106.18769  7 -2.417842  0.0462

#########################################################
model5 <- lme(daycal_CaCO3.cm2.min~condition,data=dat,random= ~1|sid,method="REML",na.action=na.omit) #change to different traits


#Diagnostics to determine if models are any good to begin with
qqPlot(residuals(model5),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal

qqnorm(model5, ~ranef(., level=1))


#check homoscedascitity by plotting residuals against fitted values
plot(model5,which=1) #want random scatter; no apparent trendline


summary(model5)
anova(model5)

Fixed effects: daycal_CaCO3.cm2.min ~ condition 
                    Value Std.Error DF   t-value p-value
(Intercept)      41.04176  4.482284  7  9.156439  0.0000
conditionbleach -11.71891  3.375153  7 -3.472111  0.0104

#########################################################
#PCA of physiology variables for GE correlation

head(dat)

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE.
#library(devtools)

 
pca <- prcomp(dat[,c(4,6,9:11)],
      center = TRUE,
      scale. = TRUE) 
summary(pca) #First 2 PCs explain 61% of the variation
str(pca)

scores=pca$x

scores[,1:2]
dat[,1:2]

library(ggfortify)
df<-dat[,c(4,6,9:11)]
pP<-autoplot(prcomp(df,center=TRUE,scale.=TRUE),data=dat,colour="condition",loadings=TRUE,loadings.colour="grey60",loadings.label=TRUE,loadings.label.size=3,loadings.label.colour="grey60")+scale_colour_manual(values=c("darkgreen","lightgreen"))+theme_bw()
pG<-autoplot(prcomp(df,center=TRUE,scale.=TRUE),data=dat,colour="sid",loadings=TRUE,loadings.colour="grey60",loadings.label=TRUE,loadings.label.size=3,loadings.label.colour="grey60")+theme_bw()

ggarrange(pP,pG,ncol=1,labels=c('a','b'),label.args=list(gp=grid::gpar(fint=4,cex=1.2)))


#1st PC explains 64% of the variation and distinguishes normal from WPS phenotypes, 2nd PC is a lot of genotype, but not totally

#########################################################

library(dabestr)

# Performing unpaired (two independent groups) analysis.
paired_mean_diff_GrossPS <- dabest(dat, condition, GrossPS_mgO2.cm2.min,
                             idx = c("healthy", "WPS"),
                             paired = TRUE,
                             id.column = sid)

# Display the results in a user-friendly format.
paired_mean_diff_GrossPS
# DABEST (Data Analysis with Bootstrap Estimation) v0.2.0
# =======================================================

# Variable: GrossPS_mgO2.cm2.min 

# Paired mean difference of bleach (n=8) minus normal (n=8)
 # -0.898 [95CI  -1.15; -0.602]


# 5000 bootstrap resamples.
# All confidence intervals are bias-corrected and accelerated.

# Produce a Cumming estimation plot.
plot(paired_mean_diff_GrossPS,rawplot.ylabel=expression(mgO["2"]/cm["2"]/min))

