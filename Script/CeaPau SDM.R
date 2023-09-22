#Species distribution models analyzing Ceanothus pauciflorus complex data
#Using GBIF data for C. pauciflorus 
#Based on code borrowed from J.Yoder
#06/24/2023
#created by Shane E. Jordan

#clear the environment
rm(list=ls())

#required libraries
library(dplyr)
library(maptools)
library(maps)
library(mapdata)
library(dismo) #dismo has the SDM analyses we'll need
library(raster)
#also -- be sure to deselect purr, tidyr, and tidyvrse from installed packages

#read and clean coordinate data from GBIF (using Ceanothus pauciflorus here)
CeaPau <- read.csv('Data/CeaPau.csv') #coords of 308 observations - note gbif() works too
CeaPau <- rename(CeaPau, longitude = ï..longitude)

#filter data for extent of interest
CeaPau <- CeaPau %>%
  filter(longitude <= -115.0) %>%
  filter(longitude >= -118.5) %>%
  filter(latitude <=34.5) %>%
  filter(latitude >= 30.5)

#create SDM functions

#make function SDMtest to check data distribution

SDMtest <- function(data) {
  data(stateMapEnv) # load the database with the U.S. state borders
  plot(c(-118.5, -115), c(30.5, 34.5), mar=par("mar"), xlab="longitude", ylab="latitude", xaxt="n", yaxt="n", type="n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="lightblue")
  map("state", xlim=c(-118.5, -115), ylim=c(30.5, 34.5), fill=T, col="cornsilk", add=T)
  points(data$longitude, data$latitude, col="darkolivegreen4", pch=20, cex=0.5)
  axis(1,las=1)
  axis(2,las=1)
  box()
  longrid = seq(-118.5,-115,0.05)
  latgrid = seq(30.5,34.5,0.05)
}

#make function SDM to create present and projected SDMs

SDM <- function (data, test) {
  x = circles(data[,c("longitude","latitude")], d=50000, lonlat=T)
  bg = spsample(x@polygons, 1000, type='random', iter=1000)
  BClim = getData("worldclim", var="bio", res=2.5, path="data/")
  dataRange = extent(-118.5,-115,30.5,34.5)
  BClim = crop(BClim, dataRange)
  writeRaster(BClim, filename="data/dataBC_2.5.grd", overwrite=T)
  BClim = brick("data/dataBC_2.5.grd")
  data_bc = extract(BClim, data$longitude, data$latitude) # for presence points
  bg_bc = extract(BClim, bg) # for pseudo-absence points
  #build data frames
  data_bc = data.frame(lon=data$longitude, lat=data$latitude, data_bc)
  bgpoints = bg@coords
  colnames(bgpoints) = c("lon","lat")
  bg_bc = data.frame(cbind(bgpoints,bg_bc))
  length(which(is.na(bg_bc$bio1))) # double-check for missing data
  bg_bc = bg_bc[!is.na(bg_bc$bio1), ] # and pull out the missing lines
  #estimate predictive model where p=presence and a=absence
  group_p = kfold(data_bc, 5) # vector of group assignments splitting the data_bc into 5 groups
  group_a = kfold(bg_bc, 5) # ditto for bg_bc
  #validation test
  test = test
  train_p = data_bc[group_p!=test, c("lon","lat")]
  train_a = bg_bc[group_a!=test, c("lon","lat")]
  test_p = data_bc[group_p==test, c("lon","lat")]
  test_a = bg_bc[group_a==test, c("lon","lat")]
  #estimate maxent model
  me = maxent(BClim, p=train_p, a=train_a)
  e =  evaluate(test_p, test_a, me, BClim) #me=model, BClim=predictor variables
  #generate the predictions
  pred_me = predict(me, BClim) 
  #make plot
  map(xlim=c(-118.5, -115), ylim=c(30.5, 34.5), fill=F, col="black", add=F)
  plot1 <- plot(pred_me, main="Maxent, raw values - Present")
  points(data$longitude, data$latitude, pch=20, cex=0.5, col="darkgreen")
  map(xlim=c(-118.5, -115), ylim=c(30.5, 34.5), fill=F, col="black", add=T)
  ProjClim <- getData("CMIP5", var="bio", res=2.5, rcp=45, model='MI', year = 70, path="data/")
  names(ProjClim) <- names(BClim)
  pfut <- predict(ProjClim, me, ext=dataRange, progress='')
  plot2 <- plot(pfut, main="Maxent, raw values - Future")
  points(data$longitude, data$latitude, pch=20, cex=0.5, col="darkgreen")
  map(xlim=c(-118.5, -115), ylim=c(30.5, 34.5), fill=F, col="black", add=T)
  pme = maxent(ProjClim, p=train_p, a=train_a)
  pe = evaluate(test_p, test_a, pme, ProjClim) #me=model, ProjClim=predictor variables
  projdata <- rasterToPoints(pfut)
  projdata <- as.data.frame(projdata)
  predata <- rasterToPoints(pred_me)
  predata <- as.data.frame(predata)
  return(list(e,me, pe, pme, predata, projdata))
  
}

#make function to extract suitability data

SDMdata <- function(data, predmin, projmin) {
  predata <- data[[5]]
  predsuitable <- predata %>% filter(layer >  predmin)
  str(predsuitable)
  projdata <- data[[6]]
  projsuitable <- projdata %>% filter(layer > projmin)
  str(projsuitable)
}

#load data
CeaPau <- read.csv('Data/CeaPau_sp.csv')

#produce SDMs
CPau1 <- SDM(CeaPau, 1)
#auc, cor, thresholds
CPau1[[1]] #AUC=0.85
CPau1[[3]] #AUC=0.86
#look at at variables with highest contribution
CPau1[[2]] #Bio11, Bio17, Bio7, Bio12
CPau1[[4]] #Bio11, Bio13, Bio7, Bio17

CPau2 <- SDM(CeaPau, 2)
#auc, cor, thresholds
CPau2[[1]] #AUC=0.82
CPau2[[3]] #AUC=0.81
#look at at variables with highest contribution
CPau2[[2]] #Bio11, Bio14, Bio12, Bio17
CPau2[[4]] #Bio11, Bio14, Bio16, Bio7

CPau3 <- SDM(CeaPau, 3)
#auc, cor, thresholds
CPau3[[1]] #AUC=0.88
CPau3[[3]] #AUC=0.88
#look at at variables with highest contribution
CPau3[[2]] #Bio11, Bio13, Bio4, Bio16
CPau3[[4]] #Bio11, Bio16, Bio14, Bio13

CPau4 <- SDM(CeaPau, 4)
#auc, cor, thresholds
CPau4[[1]] #AUC=0.87
CPau4[[3]] #AUC=0.86
#look at at variables with highest contribution
CPau4[[2]] #Bio11, Bio16, Bio7, Bio4 
CPau4[[4]] #Bio11, Bio13, Bio7, Bio4

CPau5 <- SDM(CeaPau, 5)
#auc, cor, thresholds
CPau5[[1]] #AUC=0.85
CPau5[[3]] #AUC=0.84
#look at at variables with highest contribution
CPau5[[2]] #Bio 11, Bio5, Bio16, Bio7
CPau5[[4]] #Bio11, Bio13, Bio5, Bio7

#looks like k=3 training groups fits the best model
#next look at numbers of grids above and below habitat suitability thresholds from above
SDMdata(CPau,  0.4117009 , 0.3436315)
#present=1010
#future=155

