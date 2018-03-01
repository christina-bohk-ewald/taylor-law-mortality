#######################################################################
#######################################################################
## 'taylor-law-mortality' is a program that provides the R source code
## used for the calculations in the paper: 
## Cohen, Bohk-Ewald, Rau (2018). Gompertz, Makeham, and Siler models explain 
## Taylor's law in human mortality data. Demographic Research 38(29): 773-842. 
## (c) Copyright 2018, Christina Bohk-Ewald

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program (see LICENSE.txt in the root of this source code).
## If not, see <http://www.gnu.org/licenses/>.
#######################################################################
#######################################################################

#######################################################################
#######################################################################
## Determine (relative) working paths:

## Please set 'the.script.path' to the respective github path on your device.

## Main path of 'taylor-law-mortality' on your device:
the.script.path <- c(".") 

## This (relative) path contains dumped data / results.
## Warning: this (relative) path is empty due to copyrights.
## You will need to put in mortality data yourself in order to use the R source code of 'taylor-law-mortality'. 
## Details on the required data objects are given below (see section 1). 
the.dump.path <- c("./dumps/") 

## This (relative) path contains figures / plots:
the.plot.path <- c("./figs/") 

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##1. Load required functions and data:
##Note that this R source code only works after you have run "step-1-mortality-parameter-estimation.R".

setwd(the.script.path)

source("fctn-TLmoments.R")

setwd(the.dump.path)

source("deaths_wom.R")
source("deaths_men.R")
source("exposures_wom.R")
source("exposures_men.R")

source("siler_par_wom.R")
source("siler_est_wom.R")
source("siler_par_men.R")
source("siler_est_men.R")

source("gompertz_par_wom.R")
source("gompertz_est_wom.R")
source("gompertz_par_men.R")
source("gompertz_est_men.R")

source("gompertz_est_betaUpConst_wom.R")
source("gompertz_est_betaUpConst_men.R")
source("gompertz_est_betaDownConst_wom.R")
source("gompertz_est_betaDownConst_men.R")

source("gompertz_makeham_par_wom.R")
source("gompertz_makeham_est_wom.R")
source("gompertz_makeham_par_men.R")
source("gompertz_makeham_est_men.R")

#######################################################################
#######################################################################

#######################################################################
#######################################################################

##2. Calculate TL for selected countries:

countrySelect <- c("DNK","FRA","GDR","FRG","HUN","ITA","JPN","POL","RUS","SWE","GBR","USA")
countrySelect2 <- c("Denmark","France","East Germany","West Germany","Hungary","Italy","Japan","Poland","Russia","Sweden","United Kingdom","USA")

######################
##2a: TL in observed data:
######################

mx.obs.wom <- list()
mx.obs.men <- list()

for(i in 1:length(countrySelect)){
	mx.obs.wom[[countrySelect[i]]] <- deaths.wom[[countrySelect[i]]]/exposures.wom[[countrySelect[i]]]
	mx.obs.men[[countrySelect[i]]] <- deaths.men[[countrySelect[i]]]/exposures.men[[countrySelect[i]]]
}

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.obs.wom <- list()
TL.moments.obs.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.obs.wom[[countrySelect[i]]] <- TL.moments(mx.obs.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.obs.men[[countrySelect[i]]] <- TL.moments(mx.obs.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.obs.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.obs.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.obs.wom[[countrySelect[i]]]$var.mx[which(TL.moments.obs.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.obs.men[[countrySelect[i]]]$var.mx[which(TL.moments.obs.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.obs.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.obs.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.obs.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.obs.wom[i,3] <- summary(lm(log(TL.moments.obs.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.obs.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.obs.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.obs.men[i,1] <- help.coef.men[1]
	TL.cross.age.obs.men[i,2] <- help.coef.men[2]
	TL.cross.age.obs.men[i,3] <- summary(lm(log(TL.moments.obs.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.obs.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

##############
##2b: TL in Siler:
##############

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.siler.wom <- list()
TL.moments.siler.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.siler.wom[[countrySelect[i]]] <- TL.moments(siler.est.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.siler.men[[countrySelect[i]]] <- TL.moments(siler.est.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.siler.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.siler.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.siler.wom[[countrySelect[i]]]$var.mx[which(TL.moments.siler.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.siler.men[[countrySelect[i]]]$var.mx[which(TL.moments.siler.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.siler.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.siler.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.siler.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.siler.wom[i,3] <- summary(lm(log(TL.moments.siler.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.siler.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.siler.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.siler.men[i,1] <- help.coef.men[1]
	TL.cross.age.siler.men[i,2] <- help.coef.men[2]
	TL.cross.age.siler.men[i,3] <- summary(lm(log(TL.moments.siler.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.siler.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

#################
##2c: TL in Gompertz:
#################

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.gompertz.wom <- list()
TL.moments.gompertz.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.gompertz.wom[[countrySelect[i]]] <- TL.moments(gompertz.est.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.gompertz.men[[countrySelect[i]]] <- TL.moments(gompertz.est.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.gompertz.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.gompertz.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.gompertz.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.gompertz.wom[i,3] <- summary(lm(log(TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.men[i,1] <- help.coef.men[1]
	TL.cross.age.gompertz.men[i,2] <- help.coef.men[2]
	TL.cross.age.gompertz.men[i,3] <- summary(lm(log(TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

#########################
##2d: TL in Gompertz-Makeham:
#########################

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.gompertz.makeham.wom <- list()
TL.moments.gompertz.makeham.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.gompertz.makeham.wom[[countrySelect[i]]] <- TL.moments(gompertz.makeham.est.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.gompertz.makeham.men[[countrySelect[i]]] <- TL.moments(gompertz.makeham.est.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.gompertz.makeham.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.gompertz.makeham.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.makeham.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.gompertz.makeham.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.gompertz.makeham.wom[i,3] <- summary(lm(log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.makeham.men[i,1] <- help.coef.men[1]
	TL.cross.age.gompertz.makeham.men[i,2] <- help.coef.men[2]
	TL.cross.age.gompertz.makeham.men[i,3] <- summary(lm(log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

###############################
##2e: TL in Gompertz - betaUpConst:
###############################

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.gompertz.betaUpConst.wom <- list()
TL.moments.gompertz.betaUpConst.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]] <- TL.moments(gompertz.est.betaUpConst.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]] <- TL.moments(gompertz.est.betaUpConst.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.gompertz.betaUpConst.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.gompertz.betaUpConst.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.betaUpConst.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.gompertz.betaUpConst.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.gompertz.betaUpConst.wom[i,3] <- summary(lm(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.betaUpConst.men[i,1] <- help.coef.men[1]
	TL.cross.age.gompertz.betaUpConst.men[i,2] <- help.coef.men[2]
	TL.cross.age.gompertz.betaUpConst.men[i,3] <- summary(lm(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

#################################
##2f: TL in Gompertz - betaDownConst:
#################################

age.rows <- 1:101
time.col <- as.character(1960:2009) 

TL.moments.gompertz.betaDownConst.wom <- list()
TL.moments.gompertz.betaDownConst.men <- list()

for(i in 1:length(countrySelect)){
	TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]] <- TL.moments(gompertz.est.betaDownConst.wom[[countrySelect[i]]][age.rows,time.col])
	TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]] <- TL.moments(gompertz.est.betaDownConst.men[[countrySelect[i]]][age.rows,time.col])
}

TL.cross.age.gompertz.betaDownConst.wom <- matrix(0,nrow=length(countrySelect),ncol=3)
TL.cross.age.gompertz.betaDownConst.men <- matrix(0,nrow=length(countrySelect),ncol=3)

for(i in 1:length(countrySelect)){

	print(i)

	TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx==0)] <- NA
	TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx[which(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx==0)] <- NA
	
	help.coef.wom <- coef(lm(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.betaDownConst.wom[i,1] <- help.coef.wom[1]
	TL.cross.age.gompertz.betaDownConst.wom[i,2] <- help.coef.wom[2]
	TL.cross.age.gompertz.betaDownConst.wom[i,3] <- summary(lm(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared

	help.coef.men <- coef(lm(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10)),na.action=na.omit))
	TL.cross.age.gompertz.betaDownConst.men[i,1] <- help.coef.men[1]
	TL.cross.age.gompertz.betaDownConst.men[i,2] <- help.coef.men[2]
	TL.cross.age.gompertz.betaDownConst.men[i,3] <- summary(lm(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx[1:101],base=10)~(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10))))$r.squared
}

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##Create Figures 1 through 17 of Cohen, Bohk-Ewald, Rau (2018). 

##
##Figure 1 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("MortFit-models-wom.pdf",width=15,height=15,pointsize=20,family="Helvetica")

par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))

rbPal <- colorRampPalette(c('gray','black'))
int.col <- rbPal(20)[as.numeric(cut(c(1:50),breaks=20))]

for(country in 1:12){
	plot(x=-40,y=-40,xlim=c(0,100),ylim=c(log(0.00001,base=10),log(1,base=10)),xlab="Age",ylab="",axes=FALSE,main=countrySelect2[country])
	axis(side=1,at=seq(0,100,5),labels=FALSE,lwd=1,pos=log(0.00001,base=10))
	axis(side=1,at=seq(0,100,10),labels=TRUE,lwd=3,pos=log(0.00001,base=10))
	axis(side=2,at=c(log(0.00001,base=10),log(0.0001,base=10),log(0.001,base=10),log(0.01,base=10),log(0.1,base=10),log(1,base=10)),labels=c(0.00001,0.0001,0.001,0.01,0.1,1),lwd=3,pos=0)
	for(i in 1:50){
		lines(x=0:100,y=log(mx.obs.wom[[country]][1:101,as.character(1960:2009)[i]],base=10),col=int.col[i-10])
	}
	lines(x=0:100,y=log(siler.est.wom[[country]][1:101,as.character(1960:2009)[i]],base=10),col="red",lwd=2)
	lines(x=0:100,y=log(gompertz.est.wom[[country]][1:101,as.character(1960:2009)[i]],base=10),col="green",lwd=2)
	lines(x=0:100,y=log(gompertz.makeham.est.wom[[country]][1:101,as.character(1960:2009)[i]],base=10),col="blue",lwd=2)
	legend(x=-2,y=log(1.1,base=10),c("Women","Observed, 1960-2009","Gompertz, 2009","Makeham, 2009","Siler, 2009"),
		col=c(NA,"black","green","blue","red"),bty="n",cex=0.65,lty=c(NA,1,1,1,1),lwd=2)
}

dev.off()

##
##Figure 2 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("MortFit-models-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))

rbPal <- colorRampPalette(c('gray','black'))
int.col <- rbPal(20)[as.numeric(cut(c(1:50),breaks=20))]

for(country in 1:12){
	plot(x=-40,y=-40,xlim=c(0,100),ylim=c(log(0.00001,base=10),log(1,base=10)),xlab="Age",ylab="",axes=FALSE,main=countrySelect2[country])
	axis(side=1,at=seq(0,100,5),labels=FALSE,lwd=1,pos=log(0.00001,base=10))
	axis(side=1,at=seq(0,100,10),labels=TRUE,lwd=3,pos=log(0.00001,base=10))
	axis(side=2,at=c(log(0.00001,base=10),log(0.0001,base=10),log(0.001,base=10),log(0.01,base=10),log(0.1,base=10),log(1,base=10)),labels=c(0.00001,0.0001,0.001,0.01,0.1,1),lwd=3,pos=0)
	for(i in 1:50){
		lines(x=0:100,y=log(mx.obs.men[[country]][1:101,as.character(1960:2009)[i]],base=10),col=int.col[i-10])
	}
	lines(x=0:100,y=log(siler.est.men[[country]][1:101,as.character(1960:2009)[i]],base=10),col="red",lwd=2)
	lines(x=0:100,y=log(gompertz.est.men[[country]][1:101,as.character(1960:2009)[i]],base=10),col="green",lwd=2)
	lines(x=0:100,y=log(gompertz.makeham.est.men[[country]][1:101,as.character(1960:2009)[i]],base=10),col="blue",lwd=2)
	legend(x=-2,y=log(1.1,base=10),c("Men","Observed, 1960-2009","Gompertz, 2009","Makeham, 2009","Siler, 2009"),
		col=c(NA,"black","green","blue","red"),bty="n",cex=0.65,lty=c(NA,1,1,1,1),lwd=2)

}

dev.off()

##
##Figure 3 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-observed-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.obs.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.obs.wom[i,1] + TL.cross.age.obs.wom[i,2]*log(TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.obs.wom[i,1],2), "+ ", round(TL.cross.age.obs.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.obs.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)

}

dev.off()

##
##Figure 4 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-observed-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.obs.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.obs.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.obs.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.obs.men[i,1] + TL.cross.age.obs.men[i,2]*log(TL.moments.obs.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.obs.men[i,1],2), "+ ", round(TL.cross.age.obs.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.obs.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()


##
##Figure 5 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.wom[i,1] + TL.cross.age.gompertz.wom[i,2]*log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.wom[i,1],2), "+ ", round(TL.cross.age.gompertz.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 6 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.men[i,1] + TL.cross.age.gompertz.men[i,2]*log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.men[i,1],2), "+ ", round(TL.cross.age.gompertz.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 7 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-makeham-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.makeham.wom[i,1] + TL.cross.age.gompertz.makeham.wom[i,2]*log(TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.makeham.wom[i,1],2), "+ ", round(TL.cross.age.gompertz.makeham.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.makeham.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 8 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-makeham-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.makeham.men[i,1] + TL.cross.age.gompertz.makeham.men[i,2]*log(TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.makeham.men[i,1],2), "+ ", round(TL.cross.age.gompertz.makeham.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.makeham.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()


##
##Figure 9 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-siler-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.siler.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.siler.wom[i,1] + TL.cross.age.siler.wom[i,2]*log(TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.siler.wom[i,1],2), "+ ", round(TL.cross.age.siler.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.siler.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 10 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-siler-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.siler.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.siler.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.siler.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.siler.men[i,1] + TL.cross.age.siler.men[i,2]*log(TL.moments.siler.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.siler.men[i,1],2), "+ ", round(TL.cross.age.siler.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.siler.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()


##
##Figure 11 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-scatter-sex-slope.pdf",width=15, height=15,pointsize=20,family="Helvetica")

par(mfrow=c(1,1),las=1)

plot(x=-40,y=-40,xlim=c(1.2,2.1),ylim=c(1.2,2.1),xlab="TL, slope, women",ylab="TL, slope, men",main="",axes=FALSE)
axis(side=1,at=seq(1.2,2.1,0.1),labels=TRUE,lwd=3,pos=1.2)
axis(side=2,at=seq(1.2,2.1,0.1),labels=TRUE,lwd=3,pos=1.2)
for(i in 1:12){
	text(TL.cross.age.obs.wom[i,2],TL.cross.age.obs.men[i,2],countrySelect[i],col="black",font=2)
	text(TL.cross.age.siler.wom[i,2],TL.cross.age.siler.men[i,2],countrySelect[i],col="red",font=2)
	text(TL.cross.age.gompertz.wom[i,2],TL.cross.age.gompertz.men[i,2],countrySelect[i],col="green",font=2)
	text(TL.cross.age.gompertz.makeham.wom[i,2],TL.cross.age.gompertz.makeham.men[i,2],countrySelect[i],col="blue",font=2)
}
lines(x=seq(1.2,2.1,0.1),y=seq(1.2,2.1,0.1))
legend(x=1.25,y=2.1,c("Observed","Gompertz","Makeham","Siler"),
	col=c("black","green","blue","red"),bty="n",pch=16)

dev.off()

##
##Figure 12 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-betaUpConst-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.betaUpConst.wom[i,1] + TL.cross.age.gompertz.betaUpConst.wom[i,2]*log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.betaUpConst.wom[i,1],2), "+ ", round(TL.cross.age.gompertz.betaUpConst.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.betaUpConst.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 13 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-betaUpConst-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.betaUpConst.men[i,1] + TL.cross.age.gompertz.betaUpConst.men[i,2]*log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.betaUpConst.men[i,1],2), "+ ", round(TL.cross.age.gompertz.betaUpConst.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.betaUpConst.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 14 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-betaDownConst-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.betaDownConst.wom[i,1] + TL.cross.age.gompertz.betaDownConst.wom[i,2]*log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.betaDownConst.wom[i,1],2), "+ ", round(TL.cross.age.gompertz.betaDownConst.wom[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.betaDownConst.wom[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 15 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-betaDownConst-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('yellow','red','magenta','turquoise','blue','green'))
int.col <- rbPal(20)[as.numeric(cut(c(1:101),breaks=20))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	plot(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx[1:101],base=10),col=int.col,main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=2)
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	lines(x=log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.betaDownConst.men[i,1] + TL.cross.age.gompertz.betaDownConst.men[i,2]*log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10))
	text(x=-2.3,y=-1.6,paste("log Var =", round(TL.cross.age.gompertz.betaDownConst.men[i,1],2), "+ ", round(TL.cross.age.gompertz.betaDownConst.men[i,2],2), "* log E" ))
	lines <- c(expression(paste("r"^2)), paste("       =", round(TL.cross.age.gompertz.betaDownConst.men[i,3],6)) )
	text(x=c(-2.6,-2),y=-2.4, labels=lines)
	points(x=rep(-0.4,101),y=seq(-7.2,-5.8,length=101),col=int.col)
	text(x=-0.4,y=-7.4,"Age 0",pos=1)
	text(x=-0.4,y=-5.0,"Age 100",pos=1)
}

dev.off()

##
##Figure 16 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-allScenarios-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	##betaUpConst:
	plot(log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaUpConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col="forestgreen",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3,type="l")
	##betaDownConst:
	lines(log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaDownConst.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col="tomato",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3)
	##justBeta:
	lines(log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[1:101],base=10),col="steelblue",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3)
	lines(x=log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.wom[i,1] + TL.cross.age.gompertz.wom[i,2]*log(TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[1:101],base=10),lwd=3,col="navy")
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	legend(x=-4,y=-1.8,col=c("steelblue","navy","forestgreen","tomato"),
		c("Gompertz","Gompertz, TL",expression(beta['t,up']==paste(constant)),expression(beta['t,down']==paste(constant))),
		lwd=3,cex=0.9,bty="n")
}

dev.off()

##
##Figure 17 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("TL-crossAge-gompertz-allScenarios-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
	##betaUpConst:
	plot(log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaUpConst.men[[countrySelect[i]]]$var.mx[1:101],base=10),col="forestgreen",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3,type="l")
	##betaDownConst:
	lines(log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.betaDownConst.men[[countrySelect[i]]]$var.mx[1:101],base=10),col="tomato",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3)
	##justBeta:
	lines(log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10),log(TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[1:101],base=10),col="steelblue",main=countrySelect2[i],xlim=c(-4,0),ylim=c(-8.0,-1.5),xlab="E(m)",ylab="Var(m)",axes=FALSE,lwd=3)
	lines(x=log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10),y=TL.cross.age.gompertz.men[i,1] + TL.cross.age.gompertz.men[i,2]*log(TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[1:101],base=10),lwd=3,col="navy")
	axis(side=1,at=seq(-4,0,by=1),labels=10^seq(-4,0,by=1),pos=-8)
	axis(side=2,at=seq(-8,-2,by=2),labels=10^seq(-8,-2,by=2),pos=-4)
	legend(x=-4,y=-1.8,col=c("steelblue","navy","forestgreen","tomato"),
		c("Gompertz","Gompertz, TL",expression(beta['t,up']==paste(constant)),expression(beta['t,down']==paste(constant))),
		lwd=3,cex=0.9,bty="n")
}

dev.off()
