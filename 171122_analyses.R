###################### DIVGRASS DIVERSITY-STABILITY RELATIONSHIPS / ANALYSES
###################### By Lucie Mahaut, october/Nomvember 2022
###################### Analyse relationships between EP at both the french and divgrass module scales
###################### Disentangle the role of abiotic and biotic drivers 

#### Download general libraries
library(tidyr)
library(ggplot2)
library(plyr)
library(dplyr)
library(MASS)
library(car)
library(vegan)
library(viridis)
library(rsample)
library(ggExtra)
library(data.table)
library(FD)# functional diversity
library(Weighted.Desc.Stat) # weighted variance and co
library(randomForest)
library(interactions)
library(jtools)
library(ggiraphExtra)

library(DHARMa)
library(mgcv)
library(nlme)
library(piecewiseSEM)

rm(list=ls())

#### functio to check residuals
lm.check<-function (x) {
  par(mfrow=c(2,2))
  hist(residuals(x))
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(residuals(x)~predict(x))
}

################## Download data ################## 
setwd("~/R/ASSET/Divgrass")

# environment and floristic data
dat<-readRDS("081022_data_f.rds")

################## Process the data
datX<-dat

#### divgrass pixels
datX$divgrass<-NA
datX$divgrass[is.na(datX$MODULE)]<-0
datX$divgrass[is.na(datX$divgrass)]<-1

#### select divgrass pixels
data<-droplevels(datX[datX$divgrass==1,])
dim(data)# 1244 pixels

#### remove NA in trait
data<-data[!is.na(data$flo_m),]
data<-data[!is.na(data$ sla_v),]
data<-data[!is.na(data$ la_v),]
dim(data)# 1193 pixels

#### Change covariate if necessary

data$flo_v<-sqrt(data$flo_v)
data$la_m<-log(data$la_m)
data$sm_m<-log(data$sm_m)
data$height_v<-sqrt(data$height_v)
data$la_v<-sqrt(data$la_v)
data$sm_v<-log(data$sm_v)
data$ldmc_v<-sqrt(data$ldmc_v)

dataX<-data

#### Scale covariates
data$n_tot<-(data$n_tot-mean(data$n_tot))/sd(data$n_tot)

data$GSL<-(data$GSL-mean(data$GSL))/sd(data$GSL)

data$alpha_ediv<-(data$alpha_ediv-mean(data$alpha_ediv))/sd(data$alpha_ediv)
data$alpha_rich<-(data$alpha_rich-mean(data$alpha_rich))/sd(data$alpha_rich)

data$sla_m<-(data$sla_m - mean(data$sla_m))/sd(data$sla_m)
data$height_m<-(data$height_m - mean(data$height_m))/sd(data$height_m)
data$la_m<-(data$la_m - mean(data$la_m))/sd(data$la_m)
data$sm_m<-(data$sm_m - mean(data$sm_m))/sd(data$sm_m)
data$ldmc_m<-(data$ldmc_m - mean(data$ldmc_m))/sd(data$ldmc_m)
data$lnc_m<-(data$lnc_m - mean(data$lnc_m))/sd(data$lnc_m)
data$flo_m<-(data$flo_m - mean(data$flo_m))/sd(data$flo_m)
data$sla_v<-(data$sla_v - mean(data$sla_v))/sd(data$sla_v)
data$height_v<-(data$height_v - mean(data$height_v))/sd(data$height_v)
data$la_v<-(data$la_v - mean(data$la_v))/sd(data$la_v)
data$sm_v<-(data$sm_v - mean(data$sm_v))/sd(data$sm_v)
data$ldmc_v<-(data$ldmc_v - mean(data$ldmc_v))/sd(data$ldmc_v)
data$lnc_v<-(data$lnc_v - mean(data$lnc_v))/sd(data$lnc_v)
data$flo_v<-(data$flo_v - mean(data$flo_v))/sd(data$flo_v)

data$fdis<-(data$fdis - mean(data$fdis))/sd(data$fdis)
################## Process data by habitat ################## 

################## Split by cluster
#### split (create a list)

mod.data<-split(dataX,dataX$MODULE)
mod.data.scl<-list(NA)

#### scale covariates for analysis
for (i in 1:length(mod.data))
{
  # select data
  mod.data.scl[[i]]<-mod.data[[i]]
  
  # replace NA
  mod.data.scl[[i]]$sla_v[is.na(mod.data.scl[[i]]$sla_v)]<-0
  
  mod.data.scl[[i]]$la_v<-na.roughfix(mod.data.scl[[i]]$la_v)
  mod.data.scl[[i]]$la_m<-na.roughfix(mod.data.scl[[i]]$la_m)
  
  mod.data.scl[[i]]$sm_v<-na.roughfix(mod.data.scl[[i]]$sm_v)
  mod.data.scl[[i]]$sm_m<-na.roughfix(mod.data.scl[[i]]$sm_m)
  
  mod.data.scl[[i]]$ldmc_v<-na.roughfix(mod.data.scl[[i]]$ldmc_v)
  mod.data.scl[[i]]$ldmc_m<-na.roughfix(mod.data.scl[[i]]$ldmc_m)
  
  mod.data.scl[[i]]$flo_v<-na.roughfix(mod.data.scl[[i]]$flo_v)
  mod.data.scl[[i]]$flo_m<-na.roughfix(mod.data.scl[[i]]$flo_m)
  
  mod.data.scl[[i]]$n_tot<-(mod.data.scl[[i]]$n_tot-mean(mod.data.scl[[i]]$n_tot))/sd(mod.data.scl[[i]]$n_tot)
  mod.data.scl[[i]]$GSL<-(mod.data.scl[[i]]$GSL-mean(mod.data.scl[[i]]$GSL))/sd(mod.data.scl[[i]]$GSL)
  
  mod.data.scl[[i]]$alpha_ediv<-(mod.data.scl[[i]]$alpha_ediv-mean(mod.data.scl[[i]]$alpha_ediv))/sd(mod.data.scl[[i]]$alpha_ediv)
  mod.data.scl[[i]]$alpha_rich<-(mod.data.scl[[i]]$alpha_rich-mean(mod.data.scl[[i]]$alpha_rich))/sd(mod.data.scl[[i]]$alpha_rich)
  
  mod.data.scl[[i]]$fdis<-(mod.data.scl[[i]]$fdis-mean(mod.data.scl[[i]]$fdis))/sd(mod.data.scl[[i]]$fdis)
  
    mod.data.scl[[i]]$sla_m<-(mod.data.scl[[i]]$sla_m - mean(mod.data.scl[[i]]$sla_m))/sd(mod.data.scl[[i]]$sla_m)
  mod.data.scl[[i]]$height_m<-(mod.data.scl[[i]]$height_m - mean(mod.data.scl[[i]]$height_m))/sd(mod.data.scl[[i]]$height_m)
  mod.data.scl[[i]]$la_m<-(mod.data.scl[[i]]$la_m - mean(mod.data.scl[[i]]$la_m))/sd(mod.data.scl[[i]]$la_m)
  mod.data.scl[[i]]$sm_m<-(mod.data.scl[[i]]$sm_m - mean(mod.data.scl[[i]]$sm_m))/sd(mod.data.scl[[i]]$sm_m)
  mod.data.scl[[i]]$ldmc_m<-(mod.data.scl[[i]]$ldmc_m - mean(mod.data.scl[[i]]$ldmc_m))/sd(mod.data.scl[[i]]$ldmc_m)
  mod.data.scl[[i]]$lnc_m<-(mod.data.scl[[i]]$lnc_m - mean(mod.data.scl[[i]]$lnc_m))/sd(mod.data.scl[[i]]$lnc_m)
  mod.data.scl[[i]]$flo_m<-(mod.data.scl[[i]]$flo_m - mean(mod.data.scl[[i]]$flo_m))/sd(mod.data.scl[[i]]$flo_m)
  mod.data.scl[[i]]$sla_v<-(mod.data.scl[[i]]$sla_v - mean(mod.data.scl[[i]]$sla_v))/sd(mod.data.scl[[i]]$sla_v)
  mod.data.scl[[i]]$height_v<-(mod.data.scl[[i]]$height_v - mean(mod.data.scl[[i]]$height_v))/sd(mod.data.scl[[i]]$height_v)
  mod.data.scl[[i]]$la_v<-(mod.data.scl[[i]]$la_v - mean(mod.data.scl[[i]]$la_v))/sd(mod.data.scl[[i]]$la_v)
  mod.data.scl[[i]]$sm_v<-(mod.data.scl[[i]]$sm_v - mean(mod.data.scl[[i]]$sm_v))/sd(mod.data.scl[[i]]$sm_v)
  mod.data.scl[[i]]$ldmc_v<-(mod.data.scl[[i]]$ldmc_v - mean(mod.data.scl[[i]]$ldmc_v))/sd(mod.data.scl[[i]]$ldmc_v)
  mod.data.scl[[i]]$lnc_v<-(mod.data.scl[[i]]$lnc_v - mean(mod.data.scl[[i]]$lnc_v))/sd(mod.data.scl[[i]]$lnc_v)
  mod.data.scl[[i]]$flo_v<-(mod.data.scl[[i]]$flo_v - mean(mod.data.scl[[i]]$flo_v))/sd(mod.data.scl[[i]]$flo_v)
}

#### One dataframe per habitat

calcareous_sc<-mod.data.scl[[1]]
mesic_sc<-mod.data.scl[[2]]
mountain_sc<-mod.data.scl[[3]]
ruderal_sc<-mod.data.scl[[4]]


################## Effects of habitat on EP and environmental gradients ##################
library(gridExtra)

#### Graphes EP relationships
g1<-ggplot(data, aes(MODULE,mean_kndvi,fill=MODULE))+geom_boxplot(alpha=0.5,show.legend = F)+
  scale_fill_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),axis.text.x=element_text(angle=90))+
  ylab("Productivity (kNDVI/day)") +xlab("")

g2<-ggplot(data, aes(MODULE,constancy_kndvi,fill=MODULE))+geom_boxplot(alpha=0.5,show.legend = F)+
  scale_fill_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),axis.text.x=element_text(angle=90))+
  ylab("Constancy") +xlab("")

g3<-ggplot(data, aes(MODULE,log(resistance_kndvi_dtrd),fill=MODULE))+geom_boxplot(alpha=0.5,show.legend = F)+
  scale_fill_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),axis.text.x=element_text(angle=90))+
  ylab("Resistance (log)") +xlab("")

grid.arrange(g1,g2,g3,ncol=3)

#### test
summary(aov(mean_kndvi~MODULE,data))
TukeyHSD(aov(mean_kndvi~MODULE,data))

summary(aov(constancy_kndvi~MODULE,data))
TukeyHSD(aov(constancy_kndvi~MODULE,data))

summary(aov(log(resistance_kndvi_dtrd)~MODULE,data))
TukeyHSD(aov(log(resistance_kndvi_dtrd)~MODULE,data))

#### Graphes environmenta/habitat
g1<-ggplot(dataX, aes(MODULE,GSL,fill=MODULE))+geom_boxplot(alpha=0.5,show.legend = F)+
  scale_fill_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),axis.text.x=element_text(angle=90))+
  ylab("GSL (day)") +xlab("")

g2<-ggplot(dataX, aes(MODULE,n_tot,fill=MODULE))+geom_boxplot(alpha=0.5,show.legend = F)+
  scale_fill_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),axis.text.x=element_text(angle=90))+
  ylab("N input") +xlab("")

grid.arrange(g1,g2,ncol=2)

################## Relationships between EP at both french and grassland types scales ##################

################## Compute correlation index 

#### correlations all habitat
c1<-cor.test(data$mean_kndvi,data$constancy_kndvi)
c2<-cor.test(data$mean_kndvi,log(data$resistance_kndvi))
c3<-cor.test(log(data$resistance_kndvi),data$constancy_kndvi)

#### Correlation by module
c4<-cor.test(data[data$MODULE=="Calcareous",]$mean_kndvi,data[data$MODULE=="Calcareous",]$constancy_kndvi)
c5<-cor.test(data[data$MODULE=="Calcareous",]$mean_kndvi,log(data[data$MODULE=="Calcareous",]$resistance_kndvi))
c6<-cor.test(log(data[data$MODULE=="Calcareous",]$resistance_kndvi),data[data$MODULE=="Calcareous",]$constancy_kndvi)

c7<-cor.test(data[data$MODULE=="Mesic",]$mean_kndvi,data[data$MODULE=="Mesic",]$constancy_kndvi)
c8<-cor.test(data[data$MODULE=="Mesic",]$mean_kndvi,log(data[data$MODULE=="Mesic",]$resistance_kndvi))
c9<-cor.test(log(data[data$MODULE=="Mesic",]$resistance_kndvi),data[data$MODULE=="Mesic",]$constancy_kndvi)

c10<-cor.test(data[data$MODULE=="Mountain",]$mean_kndvi,data[data$MODULE=="Mountain",]$constancy_kndvi)
c11<-cor.test(data[data$MODULE=="Mountain",]$mean_kndvi,log(data[data$MODULE=="Mountain",]$resistance_kndvi))
c12<-cor.test(log(data[data$MODULE=="Mountain",]$resistance_kndvi),data[data$MODULE=="Mountain",]$constancy_kndvi)

c13<-cor.test(data[data$MODULE=="Ruderal",]$mean_kndvi,data[data$MODULE=="Ruderal",]$constancy_kndvi)
c14<-cor.test(data[data$MODULE=="Ruderal",]$mean_kndvi,log(data[data$MODULE=="Ruderal",]$resistance_kndvi))
c15<-cor.test(log(data[data$MODULE=="Ruderal",]$resistance_kndvi),data[data$MODULE=="Ruderal",]$constancy_kndvi)

################## Graphes EP relationships
g1<-ggplot(data, aes(mean_kndvi,log(constancy_kndvi)))+
  geom_point(alpha=0.5,col="grey",show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12))+
  xlab("Productivity (kNDVI/year)") +ylab("Constancy")

g2<-ggplot(data, aes(mean_kndvi,constancy_kndvi,col=MODULE))+facet_grid(~MODULE,scales = "free_x")+
  geom_point(alpha=0.5,show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12),strip.text.x = element_text(size = 12))+
  xlab("Productivity (kNDVI/year)") +ylab("Constancy")

g3<-ggplot(data, aes(mean_kndvi,log(resistance_kndvi),col=MODULE))+
  geom_point(alpha=0.5,col="grey",show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12))+
  xlab("Productivity (kNDVI/year)") +ylab("Resistance (log)")

g4<-ggplot(data, aes(mean_kndvi,log(resistance_kndvi),col=MODULE))+facet_grid(~MODULE,scales = "free_x")+
  geom_point(alpha=0.5,show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12),strip.text.x = element_text(size = 12))+
  xlab("Productivity (kNDVI/year)") +ylab("Resistance (log)")

g5<-ggplot(data, aes(constancy_kndvi,log(resistance_kndvi),col=MODULE))+
  geom_point(alpha=0.5,col="grey",show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12))+
  xlab("Constancy") +ylab("Resistance (log)")

g6<-ggplot(data, aes(constancy_kndvi,log(resistance_kndvi),col=MODULE))+facet_grid(~MODULE,scales = "free_x")+
  geom_point(alpha=0.5,show.legend = F)+geom_smooth(method="lm", formula = y~x,col="black")+
  scale_colour_viridis(discrete = T, option = "viridis", direction = 1)+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12),strip.text.x = element_text(size = 12))+
  xlab("Constancy") +ylab("Resistance (log)")

g1
g2
g3
g4
g5
g6

################## Spatial auto-correlation : Moran test and GLS structure ################## 

################## Moran test ################## 
#### regression
mean_ns <-lm(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL + MODULE, data = data)

cons_ns <- lm(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL + MODULE , weights =MASK70_5K, data = data)


resi_ns <- lm(log(resistance_kndvi_dtrd)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL + MODULE, weights =MASK70_5K, data = data)

#### moran test
sims_mean   <- simulateResiduals(mean_ns)
testSpatialAutocorrelation(sims_mean, x = data$x, y = data$y)

sims_cons   <- simulateResiduals(cons_ns)
testSpatialAutocorrelation(sims_cons, x = data$x, y = data$y)

sims_res   <- simulateResiduals(resi_ns)
testSpatialAutocorrelation(sims_res, x = data$x, y = data$y)

################## Identify best corelation structure ################## 
#### variogram
mean_ns<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL,data=data)

cons_ns<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL,data=data)

resi_ns<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL,data=data)

plot(nlme:::Variogram(mean_ns, form = ~x +y, resType = "normalized"))
plot(nlme:::Variogram(cons_ns, form = ~x +y, resType = "normalized"))
plot(nlme:::Variogram(resi_ns, form = ~x +y, resType = "normalized"))

################## correlation structure for mean
gls_1<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL ,
              correlation = corSpatial(form = ~ x + y,type ="spherical", nugget =F ),
              data = data)

gls_2<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL ,
              correlation = corSpatial(form = ~ x + y,type ="exponential" , nugget =F ),
              data = data)

gls_3<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL ,
              correlation = corSpatial(form = ~ x + y,type ="gaussian" , nugget =F ),
              data = data)

gls_4<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL ,
              correlation = corSpatial(form = ~ x + y,type ="linear" , nugget =F ),
              data = data)

gls_5<-gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                n_tot + GSL ,
              correlation = corSpatial(form = ~ x + y,type ="rational" , nugget =F ),
              data = data)

AIC(gls_1)
AIC(gls_2) # lowest value --> the best correlation structure
AIC(gls_3)
AIC(gls_4)
AIC(gls_5)
summary(gls_2)

################## correlation structure for constancy
gls_cons1<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
             n_tot + GSL ,
           correlation = corSpatial(form = ~ x + y,type ="spherical", nugget =F ),
           data = data)

gls_cons2<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
             n_tot + GSL ,
           correlation = corSpatial(form = ~ x + y,type ="exponential" , nugget =F ),
           data = data)

gls_cons3<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
             n_tot + GSL ,
           correlation = corSpatial(form = ~ x + y,type ="gaussian" , nugget =F ),
           data = data)

gls_cons4<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
             n_tot + GSL ,
           correlation = corSpatial(form = ~ x + y,type ="linear" , nugget =F ),
           data = data)

gls_cons5<-gls(constancy_kndvi~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
             n_tot + GSL ,
           correlation = corSpatial(form = ~ x + y,type ="rational" , nugget =F ),
           data = data)

AIC(gls_cons1)
AIC(gls_cons2) # lowest value --> the best correlation structure
AIC(gls_cons3)
AIC(gls_cons4)
AIC(gls_cons5)
summary(gls_cons2)

################## correlation structure for resistance
gls_res1<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL ,
               correlation = corSpatial(form = ~ x + y,type ="spherical", nugget =F ),
               data = data)

gls_res2<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL ,
               correlation = corSpatial(form = ~ x + y,type ="exponential" , nugget =F ),
               data = data)

gls_res3<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL ,
               correlation = corSpatial(form = ~ x + y,type ="gaussian" , nugget =F ),
               data = data)

gls_res4<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL ,
               correlation = corSpatial(form = ~ x + y,type ="linear" , nugget =F ),
               data = data)

gls_res5<-gls(log(resistance_kndvi)~ sla_m + lnc_m + ldmc_m + height_m +flo_m + alpha_ediv + fdis +
                 n_tot + GSL ,
               correlation = corSpatial(form = ~ x + y,type ="rational" , nugget =F ),
               data = data)

AIC(gls_res1)
AIC(gls_res2) # lowest value --> the best correlation structure
AIC(gls_res3)
AIC(gls_res4)
AIC(gls_res5)
summary(gls_res2)

#### save the structure
csExp_mean <- corExp(form = ~ x + y, nugget = F)

################## SEM ##################
#### Correlation structure
csExp_mean <- corExp(form = ~ x + y, nugget = F)

################## calcareous
calcareous_sc$resistance_log<-log(calcareous_sc$resistance_kndvi_dtrd)

##### complete model
psem_cal <- psem(
  gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL , correlation = csExp_mean , data =calcareous_sc),
  
  gls(constancy_kndvi~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =calcareous_sc),
  
  gls(resistance_log~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =calcareous_sc),
  
  gls (sla_m ~ n_tot + GSL, data =calcareous_sc),
  gls (lnc_m ~ n_tot + GSL, data =calcareous_sc),
  gls (ldmc_m ~ n_tot + GSL, data =calcareous_sc),
  gls (height_m ~ n_tot + GSL, data =calcareous_sc),
  gls (flo_m ~ n_tot + GSL, data =calcareous_sc),
  
  gls (alpha_ediv ~ n_tot + GSL, data =calcareous_sc),
  gls (fdis ~ n_tot + GSL, data =calcareous_sc),
  
  resistance_log %~~%  constancy_kndvi,
  mean_kndvi %~~%  constancy_kndvi,
  resistance_log %~~%  mean_kndvi,
  
  fdis %~~% alpha_ediv,
  
  fdis %~~% lnc_m,
  fdis %~~% sla_m,
  fdis %~~% ldmc_m,
  fdis %~~% height_m,
  fdis %~~% flo_m,
  
  
  alpha_ediv %~~% lnc_m,
  alpha_ediv %~~% sla_m,
  alpha_ediv %~~% ldmc_m,
  alpha_ediv %~~% height_m,
  alpha_ediv %~~% flo_m,
  
  sla_m %~~% ldmc_m,
  sla_m %~~% lnc_m,
  sla_m %~~% height_m,
  sla_m %~~% flo_m,
  
  ldmc_m %~~% lnc_m,
  ldmc_m %~~% height_m,
  ldmc_m %~~% flo_m,
  
  lnc_m %~~% height_m,
  lnc_m %~~% flo_m,  
  
  height_m %~~% flo_m, 
  
  n_tot %~~% GSL)

summary(psem_cal)

################## mesic
mesic_sc$resistance_log<-log(mesic_sc$resistance_kndvi_dtrd)

##### Complete model
psem_mes <- psem(
  gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL , correlation = csExp_mean , data = mesic_sc),
  
  gls(constancy_kndvi~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =mesic_sc),
  
  gls(resistance_log~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =mesic_sc),
  
  gls (sla_m ~ n_tot + GSL, data =mesic_sc),
  gls (lnc_m ~ n_tot + GSL, data =mesic_sc),
  gls (ldmc_m ~ n_tot + GSL,data =mesic_sc),
  gls (height_m ~ n_tot + GSL,data =mesic_sc),
  gls (flo_m ~ n_tot + GSL,data =mesic_sc),
  
  gls (alpha_ediv ~ n_tot + GSL, data =mesic_sc),
  gls (fdis ~ n_tot + GSL, data =mesic_sc),
  
  resistance_log %~~%  constancy_kndvi,
  mean_kndvi %~~%  constancy_kndvi,
  resistance_log %~~%  mean_kndvi,
  
  fdis %~~% alpha_ediv,
  
  fdis %~~% lnc_m,
  fdis %~~% sla_m,
  fdis %~~% ldmc_m,
  fdis %~~% height_m,
  fdis %~~% flo_m,
  
  
  alpha_ediv %~~% lnc_m,
  alpha_ediv %~~% sla_m,
  alpha_ediv %~~% ldmc_m,
  alpha_ediv %~~% height_m,
  alpha_ediv %~~% flo_m,
  
  sla_m %~~% ldmc_m,
  sla_m %~~% lnc_m,
  sla_m %~~% height_m,
  sla_m %~~% flo_m,
  
  ldmc_m %~~% lnc_m,
  ldmc_m %~~% height_m,
  ldmc_m %~~% flo_m,
  
  lnc_m %~~% height_m,
  lnc_m %~~% flo_m,  
  
  height_m %~~% flo_m, 
  
  n_tot %~~% GSL)

summary(psem_mes)

################## Mountain

mountain_sc$resistance_log<-log(mountain_sc$resistance_kndvi_dtrd)

psem_mon <- psem(
  gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL , correlation = csExp_mean , data = mountain_sc),
  
  gls(constancy_kndvi~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =mountain_sc),
  
  gls(resistance_log~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =mountain_sc),
  
  gls (sla_m ~ n_tot + GSL,data =mountain_sc),
  gls (lnc_m ~ n_tot + GSL, data =mountain_sc),
  gls (ldmc_m ~ n_tot + GSL,data =mountain_sc),
  gls (height_m ~ n_tot + GSL, data =mountain_sc),
  gls (flo_m ~ n_tot + GSL, data =mountain_sc),
  
  gls (alpha_ediv ~ n_tot + GSL, data =mountain_sc),
  gls (fdis ~ n_tot + GSL, data =mountain_sc),
  
  resistance_log %~~%  constancy_kndvi,
  mean_kndvi %~~%  constancy_kndvi,
  resistance_log %~~%  mean_kndvi,
  
  fdis %~~% alpha_ediv,
  
  fdis %~~% lnc_m,
  fdis %~~% sla_m,
  fdis %~~% ldmc_m,
  fdis %~~% height_m,
  fdis %~~% flo_m,
  
  
  alpha_ediv %~~% lnc_m,
  alpha_ediv %~~% sla_m,
  alpha_ediv %~~% ldmc_m,
  alpha_ediv %~~% height_m,
  alpha_ediv %~~% flo_m,
  
  sla_m %~~% ldmc_m,
  sla_m %~~% lnc_m,
  sla_m %~~% height_m,
  sla_m %~~% flo_m,
  
  ldmc_m %~~% lnc_m,
  ldmc_m %~~% height_m,
  ldmc_m %~~% flo_m,
  
  lnc_m %~~% height_m,
  lnc_m %~~% flo_m,  
  
  height_m %~~% flo_m, 
  
  n_tot %~~% GSL)

summary(psem_mon)

################## ruderal
ruderal_sc$resistance_log<-log(ruderal_sc$resistance_kndvi_dtrd)

##### Complete model

psem_rud <- psem(
  gls(mean_kndvi~ sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL , correlation = csExp_mean , data = ruderal_sc),
  
  gls(constancy_kndvi~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =ruderal_sc),
  
  gls(resistance_log~  sla_m + lnc_m + ldmc_m + height_m + flo_m + alpha_ediv + fdis +
        n_tot + GSL, correlation = csExp_mean , data =ruderal_sc),
  
  gls (sla_m ~ n_tot + GSL, data =ruderal_sc),
  gls (lnc_m ~ n_tot + GSL, data =ruderal_sc),
  gls (ldmc_m ~ n_tot + GSL, data =ruderal_sc),
  gls (height_m ~ n_tot + GSL, data =ruderal_sc),
  gls (flo_m ~ n_tot + GSL, data =ruderal_sc),
  
  gls (alpha_ediv ~ n_tot + GSL, data =ruderal_sc),
  gls (fdis ~ n_tot + GSL, data =ruderal_sc),
  
  resistance_log %~~%  constancy_kndvi,
  mean_kndvi %~~%  constancy_kndvi,
  resistance_log %~~%  mean_kndvi,
  
  fdis %~~% alpha_ediv,
  
  fdis %~~% lnc_m,
  fdis %~~% sla_m,
  fdis %~~% ldmc_m,
  fdis %~~% height_m,
  fdis %~~% flo_m,
  
  
  alpha_ediv %~~% lnc_m,
  alpha_ediv %~~% sla_m,
  alpha_ediv %~~% ldmc_m,
  alpha_ediv %~~% height_m,
  alpha_ediv %~~% flo_m,
  
  sla_m %~~% ldmc_m,
  sla_m %~~% lnc_m,
  sla_m %~~% height_m,
  sla_m %~~% flo_m,
  
  ldmc_m %~~% lnc_m,
  ldmc_m %~~% height_m,
  ldmc_m %~~% flo_m,
  
  lnc_m %~~% height_m,
  lnc_m %~~% flo_m,  
  
  height_m %~~% flo_m, 
  
  n_tot %~~% GSL)

summary(psem_rud)

################## Figures SEM ################## 
library(ggpattern)

res<-read.table("171122_sem_coef.csv",h=T,sep=";") # these are coefficient extract from the SEM analysis
res$EP<-factor(res$EP,levels = c("Productivity","Constancy","Resistance"))
res$Variable<-factor(res$Variable,levels = c("CWM","Diversity",
                                             "GSL","N input"))

res$Pattern<-factor(res$Pattern,levels = c("indirect","direct"))

ggplot(res,aes(y=Percent,x=EP,fill=Variable,pattern = Pattern))+facet_grid(~Scale)+
  geom_bar_pattern(stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = colorRampPalette(c("#0066CC","#FFFFFF","#FF8C00"))(4)) +
  scale_pattern_manual(values = c(indirect = "stripe", direct = "none")) +
  labs(x = "", y = "Relative contribution", pattern = "Effect") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+
  theme_bw()+theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  

