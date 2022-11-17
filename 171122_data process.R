###################### DIVGRASS DIVERSITY-STABILITY RELATIONSHIPS / DATA PROCESS
###################### By Lucie Mahaut, October/November 2022

###################### Compute mean, sd and stability of annual EVI ;
###################### Calculate functional and taxonomic plant diversity
###################### Extract N tot and GSL from divgrass database
### Download general libraries  
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
library(FactoMineR)
library(factoextra)

###spatial libraries

library(maptools)
library(rgdal)
library(raster)
library(geosphere)
library(spData)
library(sf)
library(spaa)
library(ecospat)
library(ncdf4)
library(chron)
library(exactextractr)
library(RNetCDF)
library(devtools)
library(BHPMF)

rm(list=ls())

#### Function
mean2<-function(x){mean(x,na.rm=T)}
sd2<-function(x){sd(x,na.rm=T)}

resistance_fun<-function(x){abs(1/min(x - mean(x))/sd(x))}

#### mask 5km (i.e. raster with final resolution)
setwd("~/R/ASSET/Divgrass/")
mask5km<-raster("MASK_5KM_L93.tif")
image(mask5km)

###################################### MODIS Data PCH ##################################################

############################# Raster Prairie 500m (PCholer) Lambert 93 et mask 5km
setwd("~/R/ASSET/VI_PCHOLER_062021/")

prairie<-raster("STH_copernicus.tif")
image(prairie)

############################# Dowload data at 500m resolution

#### grasslands at 500M and KNDVI
MASK<- raster('STH_copernicus.tif')
kndvi<-read.csv("kNDVI_VIint02_FR.csv")

summary(kndvi)

dim(kndvi) # 93524 pixels of 500m with at least 70% of grasslands

#### mask70 at 500m
m      <- c(-Inf, 70, NA,  70, 100, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
MASK70 <- raster::reclassify(MASK,rclmat,progress="text")
image(MASK70)

##### ID of the cells with at least 70% grasslands in mask
m70<-rasterToPoints(MASK70)
ID<-cellFromXY(MASK,m70[,1:2])
kndvi$X<-ID

##### add coordinates to kndvi 
kndvi_df<-data.frame(xyFromCell(MASK,ID),kndvi)

colnames(kndvi_df)[1:3]<-c("lon","lat","ID500M")

############################# Aggregate at 5Km resolution
#### mask70 at 5km
MASK70_5K<-raster::aggregate(MASK70,fact=10,fun=sum,na.rm=T,filename='MASK70_5K.tif',overwrite=T)
image(MASK70_5K) # values = number of 500m pixels with at least 70% of grasslands

#### raster kndvi 500m
coordinates(kndvi_df)<-kndvi_df[,1:2]
gridded(kndvi_df)<-T

kndvi_raster<-stack(kndvi_df[,4:24])

#### aggregate at 5km and resample to match the extent of the mask
kndvi_raster5K<-aggregate(kndvi_raster,fact=10,fun=median,na.rm=T)
kndvi_raster5K<-resample(kndvi_raster5K,MASK70_5K,method="ngb")

kndvi_raster5K_f<-stack(kndvi_raster5K,MASK70_5K)

#### dataframe of annual kndvi values at 5K
kndvi_df_5K<-as.data.frame(rasterToPoints(kndvi_raster5K_f)) # 7580   23
summary(kndvi_df_5K)
kndvi_df_5K<-kndvi_df_5K[!is.na(kndvi_df_5K$X2000),]
kndvi_df_5K<-kndvi_df_5K[!is.na(kndvi_df_5K$MASK70_5K),]

dim(kndvi_df_5K)# 6535   24

############################# Compute mean, max, stability, sd and resistance of VI
##### Mean
mean_kndvi<-apply(kndvi_df_5K[,3:23],1,mean)

##### Max
max_kndvi<-apply(kndvi_df_5K[,3:23],1,max)

##### detrended NDVI by time
year=c(1:21)
dtrd.kndvi<-matrix(NA,nrow=nrow(kndvi_df_5K),ncol = 21)
for (i in 1:nrow(kndvi_df_5K))
{
  dtrd.kndvi[i,]<-residuals(lm(as.numeric(kndvi_df_5K[i,3:23])~ poly(year,2)))
}

##### SD NDVI
sd_kndvi<-apply(dtrd.kndvi,1,sd)

##### constancy (1/CV)
constancy_kndvi<-mean_kndvi/sd_kndvi

##### Resistance as the maximal deviation from baseline prod (adapted from White 2020)
# Compute function
resistance_fun<-function(x){abs(1/min(x - mean(x))/sd(x))}
data_fun<-function(x) {1999 + which(x - mean(x) == min(x - mean(x)))[1] } # extract date of max dev

## indices
resistance_kndvi<-apply(kndvi_df_5K[,3:23],1,resistance_fun)
resistance_kndvi_dtrd<-apply(dtrd.kndvi,1,resistance_fun)

date_resistance<-apply(kndvi_df_5K[,3:23],1,data_fun)

hist(date_resistance,breaks=c(1999:2020),main="",xlab="")

##### data frame
ecoprop<-data.frame(kndvi_df_5K[,c(1:2,24)],mean_kndvi,max_kndvi,sd_kndvi,constancy_kndvi,resistance_kndvi,resistance_kndvi_dtrd,date_resistance)

#### Save
setwd("~/R/ASSET/Divgrass/")
saveRDS(kndvi_df_5K,"kndvi_df_5K.rds")
saveRDS(ecoprop,"ecoprop.rds")

############################# mapping mean, constancy and resistance of KNDVI
library(ade4)

# get french boarder in wgs 84
france.contour<-getData('GADM', country='FRA', level=0)
plot(france.contour,axes = T)

# french border in l93
france.contour_l93 <- spTransform(france.contour, CRS("+init=epsg:2154"))
class(france.contour_l93) # spatial polygons dataframe
plot(france.contour_l93,axes = T)

# french  background
map <- ggplot() +
  geom_polygon(data = france.contour_l93,aes(x = long, y = lat, group = group), colour = "black", fill = NA) 

g1<-map + geom_raster(data = ecoprop, aes(x, y,fill = mean_kndvi), show.legend = F) +
  scale_fill_viridis(option="inferno",direction = 1)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

g2<-map + geom_raster(data = ecoprop, aes(x, y,fill = log(var_kndvi)), show.legend = F) +
  scale_fill_viridis(option="inferno",direction = 1)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

g3<-map + geom_raster(data = ecoprop, aes(x, y,fill = log(resistance_kndvi)), show.legend = F) +
  scale_fill_viridis(option="inferno",direction = 1)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

g4<-map + geom_raster(data = ecoprop, aes(x, y,fill = max_kndvi), show.legend = F) +
  scale_fill_viridis(option="inferno",direction = 1)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

library(gridExtra)
grid.arrange(g1,g4,g2,g3,ncol=2)

###################################### Divgrass Data ##################################################

############################# dowload data
setwd("~/R/ASSET/Divgrass/")

#### plot survey
mat_alpha<-readRDS("mat_alpha.rds") # site * species matrix
head(rownames(mat_alpha)) # com
head(colnames(mat_alpha)) # sp
dim(mat_alpha)# 19884 2008

#### plant trait
species.traits<-readRDS("species.traits.rds")
dim(species.traits) # 2596 12

#### Ntot and GSL
divgrass.dat<-readRDS("divgrass.dat.rds")
dim(divgrass.dat) # 19884     7

############################# Tranform floraison into julian data
##### replace empty by na

## helper function
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

## transform all columns
species.traits<-species.traits %>% mutate_each(funs(empty_as_na))
species.traits$floraison<-as.factor(species.traits$floraison)
summary(species.traits)

##### change month names
months_fr <- c("-avr", "-mai", "-juin", "-juil", "-août", "-sept", "-oct", "-nov", "-déc")
months_en <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
months_nb <- c("-04", "-05", "-06", "-07", "-08", "-09", "-10", "-11", "-12")

#### extract day and month
floraison_day <- substr(species.traits$floraison, 1, 2)
floraison_month <- substr(species.traits$floraison, 3, length(species.traits$floraison))

#### convert floraison date
for (i in 1:length(months_fr)){
  floraison_month[floraison_month==months_fr[i]] <- months_nb[i]
}

jjul <- paste(floraison_day, floraison_month, sep = "")
species.traits$floraison <- as.POSIXlt(jjul, format = "%d-%m")$yday

summary(species.traits)

############################# Keep mat_alpha species present in species.traits and vice et versa

##### Ordered species traits by species names
species.traits$SPECIES<-gsub(" ", "_", species.traits$SPECIES) # replace space by _
species.traits<-species.traits[order(species.traits$SPECIES),]
rownames(species.traits)<-species.traits$SPECIES

##### change species names in mat_alpha
colnames(mat_alpha)<-gsub(" ", "_", colnames(mat_alpha)) # replace space by _
mat_alphaX<-mat_alpha[,colnames(mat_alpha)%in%rownames(species.traits)]# keep only species with known traits
mat_alphaX<-mat_alphaX[! rowSums(mat_alphaX) ==0,]
dim(mat_alphaX)

species.traits<-species.traits[rownames(species.traits)%in%colnames(mat_alphaX),]
dim(species.traits)#2002 12
sp<-rownames(species.traits)

###################################### Taxonomic diversity ####################################################

dim(mat_alpha)
alpha_tdiv<-data.frame(com = rownames(mat_alpha), 
                       alpha_rich = specnumber(mat_alpha), 
                       alpha_ediv = exp(diversity(mat_alpha)) )
summary(alpha_tdiv)


###################################### Functional diversity ####################################################

############################# Compute CWM (& CWV  - not used in the paper)

##### Create dataframes per traits
sla<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
height<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
la<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
ldmc<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
sm<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
lnc<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)
flo<-matrix(NA,nrow=nrow(mat_alphaX),ncol= 2)

##### Loop
library(Weighted.Desc.Stat)

for (i in 1:nrow(mat_alphaX))# 19883 rows
{
  print(i)
  # vector of species abundance
  mu = mat_alphaX[i,] # For each pixel, a vector of species abundance (weights)
  mu = mu[mu>0] # present species
  mu = mu/sum(mu) # relative abdce
  
  # trait values
  xsla<-na.roughfix(species.traits$SLA[which(rownames(species.traits) %in% names(mu))])
  xheight<-na.roughfix(species.traits$Height[which(rownames(species.traits) %in% names(mu))])
  xsm<-na.roughfix(species.traits$SM[which(rownames(species.traits) %in% names(mu))])
  xldmc<-na.roughfix(species.traits$LDMC[which(rownames(species.traits) %in% names(mu))])
  xla<-na.roughfix(species.traits$LA[which(rownames(species.traits) %in% names(mu))])
  xlnc<-na.roughfix(species.traits$LNC_m[which(rownames(species.traits) %in% names(mu))])
  xflo<-na.roughfix(species.traits$floraison[which(rownames(species.traits) %in% names(mu))])
 
  sla[i,1]<- w.mean(xsla,mu)
  sla[i,2]<- w.var(scale(xsla),mu)
 
  height[i,1]<- w.mean(xheight,mu)
  height[i,2]<- w.var(xheight,mu)
  
  sm[i,1]<- w.mean(xsm,mu)
  sm[i,2]<- w.var(xsm,mu)
 
  la[i,1]<- w.mean(xla,mu)
  la[i,2]<- w.var(xla,mu)
  
  ldmc[i,1]<- w.mean(xldmc,mu)
  ldmc[i,2]<- w.var(xldmc,mu)
 
  lnc[i,1]<- w.mean(xlnc,mu)
  lnc[i,2]<- w.var(xlnc,mu)
  
  flo[i,1]<- w.mean(xflo,mu)
  flo[i,2]<- w.var(xflo,mu)
  
  rm(mu,xsla,xheight,xla,xldmc,xlnc,xflo,xsm)
}   

##### dataframe
colnames(sla)<-c("sla_m","sla_v")
colnames(height)<-c("height_m","height_v")
colnames(la)<-c("la_m","la_v")
colnames(sm)<-c("sm_m","sm_v")
colnames(ldmc)<-c("ldmc_m","ldmc_v")
colnames(lnc)<-c("lnc_m","lnc_v")
colnames(flo)<-c("flo_m","flo_v")

functional.div<-data.frame(com = rownames(mat_alphaX), sla,height,la,sm,ldmc,lnc,flo) 
summary(functional.div)
dim(functional.div)

############################# Fdis multitraits
library(FD)

##### keep traits 
species.traits2<-as.data.frame(species.traits[,-c(1,8:10,12)])
summary(species.traits2)
dim(species.traits2)

##### Replace NA by mean

species.traits2$Height<-na.roughfix(species.traits2$Height)
species.traits2$SLA<-na.roughfix(species.traits2$SLA)
species.traits2$SM<-na.roughfix(species.traits2$SM)
species.traits2$LDMC<-na.roughfix(species.traits2$LDMC)
species.traits2$LA<-na.roughfix(species.traits2$LA)
species.traits2$LNC_m<-na.roughfix(species.traits2$LNC_m)
species.traits2$floraison<-na.roughfix(species.traits2$floraison)

##### log transform 
species.traits2<-log10(species.traits2)

##### Compute a distance matrix between species
sp_dist<-dist(species.traits2)

##### Fdisp
dim(mat_alphaX)
fdis<-fdisp(sp_dist,mat_alphaX)
functional.div$fdis<-log(fdis$FDis)
summary(functional.div)

setwd("~/R/ASSET/Divgrass/")
saveRDS(functional.div,"111022_functional.div.rds")

###################################### Merge diversity data with N input, GSL and ecosystem properties data #############################

############################# Open diversity data
functional.div<-readRDS("111022_functional.div.rds")
head(functional.div)

#### Merge diversity and environmental (ntot and gsl) data
functional.div$com<-as.factor(functional.div$com)
alpha_tdiv$com<-alpha_tdiv$com

diversity_data_f<-divgrass.dat%>%
  left_join(alpha_tdiv)%>%
  left_join(functional.div)

summary(diversity_data_f)

############################# Match diversity with mask5km

#### From dataframe to spatial object
div_spa<-diversity_data_f

coordinates(div_spa)<-div_spa[,2:3]
crs(div_spa) <- "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs"
div_spa_5km <- spTransform(div_spa, CRS(projection(mask5km)))

image(prairie)
points(div_spa_5km)

#### Add the cells of mask5KM that correspond to the floristic releve
diversity_data_f$cells<-cellFromXY(mask5km,div_spa_5km[,2:3])
head(diversity_data_f)

#### Add 5K coordinates
coord_5K<-xyFromCell(mask5km,cell=diversity_data_f[,25])
diversity_data_f$x<-coord_5K[,1]
diversity_data_f$y<-coord_5K[,2]

#### one module per 5km cells (i.e the most frequent per cell)
module<-diversity_data_f[,c(4,25)]
module$MODULE<-as.factor(module$MODULE)
module$cells<-as.factor(module$cells)

dim(module)
head(unique(module))
# count number of each module in each cell
n<-module%>%
  count(cells,MODULE)

# keep the most frequent module
mfq<-matrix(NA,nrow=length(unique(n$cells)),ncol=3)

for (i in 1:length(unique(n$cells)))
{
  print(i)
  test<-droplevels(n[n$cells==unique(n$cells)[i],])
  
  mfq[i,1]<-as.numeric(as.character(test$cells[1]))
  mfq[i,2]<-as.numeric(as.character(test[test$n==max(test$n),2][1]))
  mfq[i,3]<-as.numeric(as.character(test[test$n==max(test$n),3][1]))
  
  rm(test)
}

mfq[1:15,]

mfqn<-as.data.frame(mfq)
colnames(mfqn)<-c("cells","MODULE","n")
mfqn$cells<-as.factor(mfqn$cells)

#### aggregate diversity data per 5km cell
diversity_data_5K<-stats::aggregate(diversity_data_f[,c(5:24,26:27)],
                             list(cells= diversity_data_f$cells),mean)

head(diversity_data_5K)
dim(diversity_data_5K)
summary(diversity_data_5K)
diversity_data_5K$cells<-as.factor(diversity_data_5K$cells)

#### Add module
diversity_data_5K<-diversity_data_5K%>%
  left_join(mfqn[,1:2])

diversity_data_5K$MODULE[diversity_data_5K$MODULE==1]<-"Calcareous"
diversity_data_5K$MODULE[diversity_data_5K$MODULE==2]<-"Mountain"
diversity_data_5K$MODULE[diversity_data_5K$MODULE==3]<-"Mesic"
diversity_data_5K$MODULE[diversity_data_5K$MODULE==5]<-"Ruderal"

diversity_data_5K$MODULE<-as.factor(diversity_data_5K$MODULE)

#### Map
g1<-map + geom_tile(data = diversity_data_5K, aes(x, y,fill = MODULE), show.legend = T) +
  scale_fill_viridis(discrete = T,option="viridis",direction = 1,alpha=0.7)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
g1

############################# final dataset
#### download dataset
ep<-readRDS("ecoprop.rds")

#### add cell number from mask 5KM
ep$cells<-cellFromXY(mask5km,ep[,1:2])
ep$cells<-as.factor(ep$cells)

#### join env ep and diversity data
data_f<-ep%>%
  left_join(diversity_data_5K)
summary(data_f)

dim(data_f)

#### save
saveRDS(data_f,"081022_process_data.rds")

