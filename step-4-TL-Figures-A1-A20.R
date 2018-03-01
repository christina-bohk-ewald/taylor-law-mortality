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

##1. Create dataframe: TLinMortalityModels (that is also available as csv file in Supplementary material of Cohen, Bohk-Ewald, Rau (2018)):
##Note that this R source code only works after you have run "step-2-TL-Figures-1-17.R".

country <- matrix(rep(c(1:12),times=1,each=50),nrow=600,ncol=1)
sex1 <- matrix(rep(c(1),times=12,each=50),nrow=600,ncol=1)
sex2 <- matrix(rep(c(2),times=12,each=50),nrow=600,ncol=1)
year <- matrix(rep(1960:2009,times=12),nrow=600,ncol=1)
TLa.wom <-  matrix(rep(TL.cross.age.obs.wom[,1],times=12,each=50),nrow=600,ncol=1)
TLa.men <-  matrix(rep(TL.cross.age.obs.men[,1],times=12,each=50),nrow=600,ncol=1)
TLb.wom <-  matrix(rep(TL.cross.age.obs.wom[,2],times=12,each=50),nrow=600,ncol=1)
TLb.men <-  matrix(rep(TL.cross.age.obs.men[,2],times=12,each=50),nrow=600,ncol=1)
GTLa.wom <-  matrix(rep(TL.cross.age.gompertz.wom[,1],times=12,each=50),nrow=600,ncol=1)
GTLa.men <-  matrix(rep(TL.cross.age.gompertz.men[,1],times=12,each=50),nrow=600,ncol=1)
GTLb.wom <-  matrix(rep(TL.cross.age.gompertz.wom[,2],times=12,each=50),nrow=600,ncol=1)
GTLb.men <-  matrix(rep(TL.cross.age.gompertz.men[,2],times=12,each=50),nrow=600,ncol=1)
GMTLa.wom <-  matrix(rep(TL.cross.age.gompertz.makeham.wom[,1],times=12,each=50),nrow=600,ncol=1)
GMTLa.men <-  matrix(rep(TL.cross.age.gompertz.makeham.men[,1],times=12,each=50),nrow=600,ncol=1)
GMTLb.wom <-  matrix(rep(TL.cross.age.gompertz.makeham.wom[,2],times=12,each=50),nrow=600,ncol=1)
GMTLb.men <-  matrix(rep(TL.cross.age.gompertz.makeham.men[,2],times=12,each=50),nrow=600,ncol=1)
STLa.wom <-  matrix(rep(TL.cross.age.siler.wom[,1],times=12,each=50),nrow=600,ncol=1)
STLa.men <-  matrix(rep(TL.cross.age.siler.men[,1],times=12,each=50),nrow=600,ncol=1)
STLb.wom <-  matrix(rep(TL.cross.age.siler.wom[,2],times=12,each=50),nrow=600,ncol=1)
STLb.men <-  matrix(rep(TL.cross.age.siler.men[,2],times=12,each=50),nrow=600,ncol=1)

Gbetatwom <- 0 
GMtwom <- 0
Gbetatmen <- 0 
GMtmen <- 0

GMctwom <- 0 
GMbetatwom <- 0 
GMMtwom <- 0
GMctmen <- 0 
GMbetatmen <- 0 
GMMtmen <- 0

Salphatwom <- 0
Sbeta1twom <- 0 
Sctwom <- 0 
Sbeta2twom <- 0
SMtwom <- 0
Salphatmen <- 0
Sbeta1tmen <- 0 
Sctmen <- 0 
Sbeta2tmen <- 0
SMtmen <- 0

for(i in 1:12){
        Gbetatwom <- c(Gbetatwom,gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),2])
        GMtwom <- c(GMtwom,gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),1])
        Gbetatmen <- c(Gbetatmen,gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),2])
        GMtmen <- c(GMtmen,gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),1])

        GMctwom <- c(GMctwom,gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),1])
        GMbetatwom <- c(GMbetatwom,gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),3])
        GMMtwom <- c(GMMtwom,gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),2])
        GMctmen <- c(GMctmen,gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),1])
        GMbetatmen <- c(GMbetatmen,gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),3])
        GMMtmen <- c(GMMtmen,gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),2])

        Salphatwom <- c(Salphatwom,siler.par.wom[[countrySelect[i]]][as.character(1960:2009),1])
        Sbeta1twom <- c(Sbeta1twom,siler.par.wom[[countrySelect[i]]][as.character(1960:2009),2])
        Sctwom <- c(Sctwom,siler.par.wom[[countrySelect[i]]][as.character(1960:2009),3])
        Sbeta2twom <- c(Sbeta2twom,siler.par.wom[[countrySelect[i]]][as.character(1960:2009),5])
        SMtwom <- c(SMtwom,siler.par.wom[[countrySelect[i]]][as.character(1960:2009),4])

        Salphatmen <- c(Salphatmen,siler.par.men[[countrySelect[i]]][as.character(1960:2009),1])
        Sbeta1tmen <- c(Sbeta1tmen,siler.par.men[[countrySelect[i]]][as.character(1960:2009),2])
        Sctmen <- c(Sctmen,siler.par.men[[countrySelect[i]]][as.character(1960:2009),3])
        Sbeta2tmen <- c(Sbeta2tmen,siler.par.men[[countrySelect[i]]][as.character(1960:2009),5])
        SMtmen <- c(SMtmen,siler.par.men[[countrySelect[i]]][as.character(1960:2009),4])
}

G.betat.wom <- matrix(Gbetatwom[-1],nrow=600,ncol=1)
G.Mt.wom <- matrix(GMtwom[-1],nrow=600,ncol=1)
G.betat.men <- matrix(Gbetatmen[-1],nrow=600,ncol=1)
G.Mt.men <- matrix(GMtmen[-1],nrow=600,ncol=1)

GM.ct.wom <- matrix(GMctwom[-1],nrow=600,ncol=1)
GM.betat.wom <- matrix(GMbetatwom[-1],nrow=600,ncol=1)
GM.Mt.wom <- matrix(GMMtwom[-1],nrow=600,ncol=1)
GM.ct.men <- matrix(GMctmen[-1],nrow=600,ncol=1)
GM.betat.men <- matrix(GMbetatmen[-1],nrow=600,ncol=1)
GM.Mt.men <- matrix(GMMtmen[-1],nrow=600,ncol=1)

S.alphat.wom <- matrix(Salphatwom[-1],nrow=600,ncol=1)
S.beta1t.wom <- matrix(Sbeta1twom[-1],nrow=600,ncol=1)
S.ct.wom <- matrix(Sctwom[-1],nrow=600,ncol=1)
S.beta2t.wom <- matrix(Sbeta2twom[-1],nrow=600,ncol=1)
S.Mt.wom <- matrix(SMtwom[-1],nrow=600,ncol=1)
S.alphat.men <- matrix(Salphatmen[-1],nrow=600,ncol=1)
S.beta1t.men <- matrix(Sbeta1tmen[-1],nrow=600,ncol=1)
S.ct.men <- matrix(Sctmen[-1],nrow=600,ncol=1)
S.beta2t.men <- matrix(Sbeta2tmen[-1],nrow=600,ncol=1)
S.Mt.men <- matrix(SMtmen[-1],nrow=600,ncol=1)

x <- 1:50 #1960:2009 -> due to interpretation we take 1 to 50

GbetatInterceptwom <- 0
GbetatSlopewom <- 0
GMtInterceptwom <- 0
GMtSlopewom <- 0
GbetatInterceptmen <- 0
GbetatSlopemen <- 0
GMtInterceptmen <- 0
GMtSlopemen <- 0

GMctInterceptwom <- 0
GMctSlopewom <- 0
GMbetatInterceptwom <- 0
GMbetatSlopewom <- 0
GMMtInterceptwom <- 0
GMMtSlopewom <- 0
GMctInterceptmen <- 0
GMctSlopemen <- 0
GMbetatInterceptmen <- 0
GMbetatSlopemen <- 0
GMMtInterceptmen <- 0
GMMtSlopemen <- 0

SalphatInterceptwom <- 0
SalphatSlopewom <- 0
Sbeta1tInterceptwom <- 0
Sbeta1tSlopewom <- 0
SctInterceptwom <- 0
SctSlopewom <- 0
Sbeta2tInterceptwom <- 0
Sbeta2tSlopewom <- 0
SMtInterceptwom <- 0
SMtSlopewom <- 0
SalphatInterceptmen <- 0
SalphatSlopemen <- 0
Sbeta1tInterceptmen <- 0
Sbeta1tSlopemen <- 0
SctInterceptmen <- 0
SctSlopemen <- 0
Sbeta2tInterceptmen <- 0
Sbeta2tSlopemen <- 0
SMtInterceptmen <- 0
SMtSlopemen <- 0


for(i in 1:12){
        GbetatInterceptwom[i] <- lm(gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        GbetatSlopewom[i] <- lm(gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]
        GMtInterceptwom[i] <- lm(gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        GMtSlopewom[i] <- lm(gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]
        GbetatInterceptmen[i] <- lm(gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        GbetatSlopemen[i] <- lm(gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]
        GMtInterceptmen[i] <- lm(gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        GMtSlopemen[i] <- lm(gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]

        GMctInterceptwom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        GMctSlopewom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]
        GMbetatInterceptwom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[1]
        GMbetatSlopewom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[2]
        GMMtInterceptwom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        GMMtSlopewom[i] <- lm(gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]
        GMctInterceptmen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        GMctSlopemen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]
        GMbetatInterceptmen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[1]
        GMbetatSlopemen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[2]
        GMMtInterceptmen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        GMMtSlopemen[i] <- lm(gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]

        SalphatInterceptwom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        SalphatSlopewom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]
        Sbeta1tInterceptwom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        Sbeta1tSlopewom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]
        SctInterceptwom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[1]
        SctSlopewom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[2]
        Sbeta2tInterceptwom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),5] ~  x )$coef[1]
        Sbeta2tSlopewom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),5] ~  x )$coef[2]
        SMtInterceptwom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),4] ~  x )$coef[1]
        SMtSlopewom[i] <- lm(siler.par.wom[[countrySelect[i]]][as.character(1960:2009),4] ~  x )$coef[2]
        SalphatInterceptmen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[1]
        SalphatSlopemen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),1] ~  x )$coef[2]
        Sbeta1tInterceptmen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[1]
        Sbeta1tSlopemen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),2] ~  x )$coef[2]
        SctInterceptmen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[1]
        SctSlopemen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),3] ~  x )$coef[2]
        Sbeta2tInterceptmen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),5] ~  x )$coef[1]
        Sbeta2tSlopemen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),5] ~  x )$coef[2]
        SMtInterceptmen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),4] ~  x )$coef[1]
        SMtSlopemen[i] <- lm(siler.par.men[[countrySelect[i]]][as.character(1960:2009),4] ~  x )$coef[2]
}

G.betat.Intercept.wom <- matrix(rep(GbetatInterceptwom,times=12,each=50),nrow=600,ncol=1)
G.betat.Slope.wom <- matrix(rep(GbetatSlopewom,times=12,each=50),nrow=600,ncol=1)
G.Mt.Intercept.wom <- matrix(rep(GMtInterceptwom,times=12,each=50),nrow=600,ncol=1)
G.Mt.Slope.wom <- matrix(rep(GMtSlopewom,times=12,each=50),nrow=600,ncol=1)
G.betat.Intercept.men <- matrix(rep(GbetatInterceptmen,times=12,each=50),nrow=600,ncol=1)
G.betat.Slope.men <- matrix(rep(GbetatSlopemen,times=12,each=50),nrow=600,ncol=1)
G.Mt.Intercept.men <- matrix(rep(GMtInterceptmen,times=12,each=50),nrow=600,ncol=1)
G.Mt.Slope.men <- matrix(rep(GMtSlopemen,times=12,each=50),nrow=600,ncol=1)

GM.ct.Intercept.wom <- matrix(rep(GMctInterceptwom,times=12,each=50),nrow=600,ncol=1)
GM.ct.Slope.wom <- matrix(rep(GMctSlopewom,times=12,each=50),nrow=600,ncol=1)
GM.betat.Intercept.wom <- matrix(rep(GMbetatInterceptwom,times=12,each=50),nrow=600,ncol=1)
GM.betat.Slope.wom <- matrix(rep(GMbetatSlopewom,times=12,each=50),nrow=600,ncol=1)
GM.Mt.Intercept.wom <- matrix(rep(GMMtInterceptwom,times=12,each=50),nrow=600,ncol=1)
GM.Mt.Slope.wom <- matrix(rep(GMMtSlopewom,times=12,each=50),nrow=600,ncol=1)
GM.ct.Intercept.men <- matrix(rep(GMctInterceptmen,times=12,each=50),nrow=600,ncol=1)
GM.ct.Slope.men <- matrix(rep(GMctSlopemen,times=12,each=50),nrow=600,ncol=1)
GM.betat.Intercept.men <- matrix(rep(GMbetatInterceptmen,times=12,each=50),nrow=600,ncol=1)
GM.betat.Slope.men <- matrix(rep(GMbetatSlopemen,times=12,each=50),nrow=600,ncol=1)
GM.Mt.Intercept.men <- matrix(rep(GMMtInterceptmen,times=12,each=50),nrow=600,ncol=1)
GM.Mt.Slope.men <- matrix(rep(GMMtSlopemen,times=12,each=50),nrow=600,ncol=1)

S.alphat.Intercept.wom <- matrix(rep(SalphatInterceptwom,times=12,each=50),nrow=600,ncol=1)
S.alphat.Slope.wom <- matrix(rep(SalphatSlopewom,times=12,each=50),nrow=600,ncol=1)
S.beta1t.Intercept.wom <- matrix(rep(Sbeta1tInterceptwom,times=12,each=50),nrow=600,ncol=1)
S.beta1t.Slope.wom <- matrix(rep(Sbeta1tSlopewom,times=12,each=50),nrow=600,ncol=1)
S.ct.Intercept.wom <- matrix(rep(SctInterceptwom,times=12,each=50),nrow=600,ncol=1)
S.ct.Slope.wom <- matrix(rep(SctSlopewom,times=12,each=50),nrow=600,ncol=1)
S.beta2t.Intercept.wom <- matrix(rep(Sbeta2tInterceptwom,times=12,each=50),nrow=600,ncol=1)
S.beta2t.Slope.wom <- matrix(rep(Sbeta2tSlopewom,times=12,each=50),nrow=600,ncol=1)
S.Mt.Intercept.wom <- matrix(rep(SMtInterceptwom,times=12,each=50),nrow=600,ncol=1)
S.Mt.Slope.wom <- matrix(rep(SMtSlopewom,times=12,each=50),nrow=600,ncol=1)
S.alphat.Intercept.men <- matrix(rep(SalphatInterceptmen,times=12,each=50),nrow=600,ncol=1)
S.alphat.Slope.men <- matrix(rep(SalphatSlopemen,times=12,each=50),nrow=600,ncol=1)
S.beta1t.Intercept.men <- matrix(rep(Sbeta1tInterceptmen,times=12,each=50),nrow=600,ncol=1)
S.beta1t.Slope.men <- matrix(rep(Sbeta1tSlopemen,times=12,each=50),nrow=600,ncol=1)
S.ct.Intercept.men <- matrix(rep(SctInterceptmen,times=12,each=50),nrow=600,ncol=1)
S.ct.Slope.men <- matrix(rep(SctSlopemen,times=12,each=50),nrow=600,ncol=1)
S.beta2t.Intercept.men <- matrix(rep(Sbeta2tInterceptmen,times=12,each=50),nrow=600,ncol=1)
S.beta2t.Slope.men <- matrix(rep(Sbeta2tSlopemen,times=12,each=50),nrow=600,ncol=1)
S.Mt.Intercept.men <- matrix(rep(SMtInterceptmen,times=12,each=50),nrow=600,ncol=1)
S.Mt.Slope.men <- matrix(rep(SMtSlopemen,times=12,each=50),nrow=600,ncol=1)

df.women <- cbind(country,sex1,year,TLa.wom,TLb.wom,
                G.betat.wom,G.Mt.wom,
                GTLa.wom,GTLb.wom,
                G.betat.Intercept.wom,G.betat.Slope.wom,
                G.Mt.Intercept.wom,G.Mt.Slope.wom,
                GM.ct.wom,GM.betat.wom,GM.Mt.wom,
                GMTLa.wom,GMTLb.wom,
                GM.ct.Intercept.wom,GM.ct.Slope.wom,
                GM.betat.Intercept.wom,GM.betat.Slope.wom,
                GM.Mt.Intercept.wom,GM.Mt.Slope.wom,
                S.alphat.wom,S.beta1t.wom,S.ct.wom,S.beta2t.wom,S.Mt.wom,
                STLa.wom,STLb.wom,
                S.alphat.Intercept.wom,S.alphat.Slope.wom,
                S.beta1t.Intercept.wom,S.beta1t.Slope.wom,
                S.ct.Intercept.wom,S.ct.Slope.wom,
                S.beta2t.Intercept.wom,S.beta2t.Slope.wom,
                S.Mt.Intercept.wom,S.Mt.Slope.wom)

head(df.women)

df.men <- cbind(country,sex2,year,TLa.men,TLb.men,
                G.betat.men,G.Mt.men,
                GTLa.men,GTLb.men,
                G.betat.Intercept.men,G.betat.Slope.men,
                G.Mt.Intercept.men,G.Mt.Slope.men,
                GM.ct.men,GM.betat.men,GM.Mt.men,
                GMTLa.men,GMTLb.men,
                GM.ct.Intercept.men,GM.ct.Slope.men,
                GM.betat.Intercept.men,GM.betat.Slope.men,
                GM.Mt.Intercept.men,GM.Mt.Slope.men,
                S.alphat.men,S.beta1t.men,S.ct.men,S.beta2t.men,S.Mt.men,
                STLa.men,STLb.men,
                S.alphat.Intercept.men,S.alphat.Slope.men,
                S.beta1t.Intercept.men,S.beta1t.Slope.men,
                S.ct.Intercept.men,S.ct.Slope.men,
                S.beta2t.Intercept.men,S.beta2t.Slope.men,
                S.Mt.Intercept.men,S.Mt.Slope.men)

head(df.men)

colnames(df.women) <- colnames(df.men) <- c("country","sex","year","TLa","TLb",
                        "G.betat","G.Mt",
                        "G.TLa","G.TLb",
                        "G.betat.Intercept","G.betat.Slope",
                        "G.Mt.Intercept","G.Mt.Slope",
                        "GM.ct","GM.betat","GM.Mt",
                        "GM.TLa","GM.TLb",
                        "GM.ct.Intercept","GM.ct.Slope",
                        "GM.betat.Intercept","GM.betat.Slope",
                        "GM.Mt.Intercept","GM.Mt.Slope",
                        "S.alphat","S.beta1t","S.ct","S.beta2t","S.Mt",
                        "S.TLa","S.TLb",
                        "S.alphat.Intercept","S.alphat.Slope",
                        "S.beta1t.Intercept","S.beta1t.Slope",
                        "S.ct.Intercept","S.ct.Slope",
                        "S.beta2t.Intercept","S.beta2t.Slope",
                        "S.Mt.Intercept","S.Mt.Slope")
                                  
TLinMortalityModels <-  as.data.frame(rbind(df.women,df.men))                                               

head(TLinMortalityModels)

#setwd(the.dump.path)
#write.csv(x=TLinMortalityModels,file="TLinMortalityModels.csv")

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##2. Create Figures A1 through A20 of Cohen, Bohk-Ewald, Rau (2018). 

##
##Figure A1 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("G-beta-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
        text(x=1963,y=0.15,expression(beta[t]))
        text(x=1985,y=0.15,paste(" = ", round(G.betat.Intercept.wom[i*50],5), "+ ", round(G.betat.Slope.wom[i*50],5), "* t" ))
	  lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.135, labels=lines)
}

dev.off()

##
##Figure A2 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("G-beta-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.15,expression(beta[t]))
        text(x=1985,y=0.15,paste(" = ", round(G.betat.Intercept.men[i*50],5), "+ ", round(G.betat.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.135, labels=lines)
}

dev.off()

##
##Figure A3 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("G-M-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.par.wom[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))       
	  text(x=1985,y=95,paste("  = ", round(G.Mt.Intercept.wom[i*50],5), "+ ", round(G.Mt.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()

##
##Figure A4 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("G-M-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.par.men[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))  
        text(x=1985,y=95,paste("  = ", round(G.Mt.Intercept.men[i*50],5), "+ ", round(G.Mt.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()

##
##Figure A5 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("GM-c-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.002),xlab="Calendar year",ylab=expression(c[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.002,by=0.0005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.002,by=0.001),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.002,expression(c[t]))  
        text(x=1985,y=0.002,paste("  = ", round(GM.ct.Intercept.wom[i*50],5), "+ ", round(GM.ct.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.0017, labels=lines)
}

dev.off()

##
##Figure A6 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("GM-c-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.002),xlab="Calendar year",ylab=expression(c[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.002,by=0.0005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.002,by=0.001),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.002,expression(c[t]))
        text(x=1985,y=0.002,paste("  = ", round(GM.ct.Intercept.men[i*50],5), "+ ", round(GM.ct.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.0017, labels=lines)
}

dev.off()

##
##Figure A7 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("GM-beta-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),3]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),3],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.15,expression(beta[t]))
        text(x=1985,y=0.15,paste(" = ", round(GM.betat.Intercept.wom[i*50],5), "+ ", round(GM.betat.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.135, labels=lines)
}

dev.off()

##
##Figure A8 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("GM-beta-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),3]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),3],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.15,expression(beta[t]))
        text(x=1985,y=0.15,paste(" = ", round(GM.betat.Intercept.men[i*50],5), "+ ", round(GM.betat.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.135, labels=lines)
}

dev.off()

##
##Figure A9 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("GM-M-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.wom[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))
        text(x=1985,y=95,paste("  = ", round(GM.Mt.Intercept.wom[i*50],5), "+ ", round(GM.Mt.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()

##
##Figure A10 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)
        
dev.off()

pdf("GM-M-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=gompertz.makeham.par.men[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))
        text(x=1985,y=95,paste("  = ", round(GM.Mt.Intercept.men[i*50],5), "+ ", round(GM.Mt.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()

##
##Figure A11 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)        

dev.off()

pdf("S-a-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.wom[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.wom[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.04),xlab="Calendar year",ylab=expression(alpha[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.04,by=0.005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.04,by=0.01),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.04,expression(alpha[t]))
        text(x=1985,y=0.04,paste(" = ", round(S.alphat.Intercept.wom[i*50],5), "+ ", round(S.alphat.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.036, labels=lines)
}

dev.off()

##
##Figure A12 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)      

dev.off()

pdf("S-a-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.men[[countrySelect[i]]][as.character(1960:2009),1]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.men[[countrySelect[i]]][as.character(1960:2009),1],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.04),xlab="Calendar year",ylab=expression(alpha[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.04,by=0.005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.04,by=0.01),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.04,expression(alpha[t]))
        text(x=1985,y=0.04,paste("  = ", round(S.alphat.Intercept.men[i*50],5), "+ ", round(S.alphat.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.036, labels=lines)
}

dev.off()

##
##Figure A13 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)  
        
dev.off()

pdf("S-beta1-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.wom[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.wom[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(2,6),xlab="Calendar year",ylab=expression(beta[1*t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=2)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=2)
        axis(side=2,at=seq(2,6,by=0.5),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(2,6,by=1),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=6,expression(beta[1*t]))
        text(x=1985,y=6,paste(" = ", round(S.beta1t.Intercept.wom[i*50],5), "+ ", round(S.beta1t.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=5.7, labels=lines)
}

dev.off()

##
##Figure A14 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)  

dev.off()

pdf("S-beta1-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.men[[countrySelect[i]]][as.character(1960:2009),2]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.men[[countrySelect[i]]][as.character(1960:2009),2],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(2,6),xlab="Calendar year",ylab=expression(beta[1*t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=2)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=2)
        axis(side=2,at=seq(2,6,by=0.5),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(2,6,by=1),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=6,expression(beta[1*t]))
        text(x=1985,y=6,paste("  = ", round(S.beta1t.Intercept.men[i*50],5), "+ ", round(S.beta1t.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=5.7, labels=lines)
}

dev.off()

##
##Figure A15 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)  

dev.off()

pdf("S-c-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.wom[[countrySelect[i]]][as.character(1960:2009),3]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.wom[[countrySelect[i]]][as.character(1960:2009),3],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.002),xlab="Calendar year",ylab=expression(c[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.002,by=0.0005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.002,by=0.001),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.002,expression(c[t]))
        text(x=1985,y=0.002,paste("  = ", round(S.ct.Intercept.wom[i*50],5), "+ ", round(S.ct.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.0017, labels=lines)
}

dev.off()

##
##Figure A16 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("S-c-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.men[[countrySelect[i]]][as.character(1960:2009),3]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.men[[countrySelect[i]]][as.character(1960:2009),3],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0,0.002),xlab="Calendar year",ylab=expression(c[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0)
        axis(side=2,at=seq(0,0.002,by=0.0005),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0,0.002,by=0.001),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.002,expression(c[t]))
        text(x=1985,y=0.002,paste("  = ", round(S.ct.Intercept.men[i*50],5), "+ ", round(S.ct.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.0017, labels=lines)
}

dev.off()

##
##Figure A17 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("S-beta2-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.wom[[countrySelect[i]]][as.character(1960:2009),5]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.wom[[countrySelect[i]]][as.character(1960:2009),5],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[2*t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
        text(x=1963,y=0.15,expression(beta[2*t]))
        text(x=1985,y=0.15,paste("  = ", round(S.beta2t.Intercept.wom[i*50],5), "+ ", round(S.beta2t.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.135, labels=lines)
}

dev.off()

##
##Figure A18 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("S-beta2-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.men[[countrySelect[i]]][as.character(1960:2009),5]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.men[[countrySelect[i]]][as.character(1960:2009),5],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(0.05,0.15),xlab="Calendar year",ylab=expression(beta[2*t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=0.05)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=0.05)
        axis(side=2,at=seq(0.05,0.15,by=0.01),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(0.05,0.15,by=0.05),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=0.15,expression(beta[2*t]))
        text(x=1985,y=0.15,paste("  = ", round(S.beta2t.Intercept.men[i*50],5), "+ ", round(S.beta2t.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=0.13, labels=lines)
}

dev.off()

##
##Figure A19 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("S-M-wom.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.wom[[countrySelect[i]]][as.character(1960:2009),4]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.wom[[countrySelect[i]]][as.character(1960:2009),4],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))
        text(x=1985,y=95,paste("  = ", round(S.Mt.Intercept.wom[i*50],5), "+ ", round(S.Mt.Slope.wom[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()

##
##Figure A20 in Cohen, Bohk-Ewald, Rau (2018):
##

setwd(the.plot.path)

dev.off()

pdf("S-M-men.pdf",width=15, height=15,pointsize=20,family="Helvetica")

rbPal <- colorRampPalette(c('lightgray','black'))
int.col <- rbPal(50)[as.numeric(cut(c(1:50),breaks=50))]
par(mfrow=c(3,4),las=1,oma=c(0,0,0,0),mar=c(4.1,4.1,1.1,0.1))
for(i in 1:12){
        current.y <- siler.par.men[[countrySelect[i]]][as.character(1960:2009),4]
        current.model <- lm(current.y ~ x)
        plot(x=1960:2009,y=siler.par.men[[countrySelect[i]]][as.character(1960:2009),4],col=int.col,main=countrySelect2[i],xlim=c(1960,2010),ylim=c(70,95),xlab="Calendar year",ylab=expression(M[t]),axes=FALSE,lwd=2)
        axis(side=1,at=seq(1960,2010,by=5),labels=FALSE,lwd=1,pos=70)
        axis(side=1,at=seq(1960,2010,by=10),labels=TRUE,lwd=3,pos=70)
        axis(side=2,at=seq(70,95,by=1),labels=FALSE,lwd=1,pos=1960)
        axis(side=2,at=seq(70,95,by=5),labels=TRUE,lwd=3,pos=1960)
        lines(x=1960:2009,y=current.model$"fitted.values")
	  text(x=1963,y=95,expression(M[t]))
        text(x=1985,y=95,paste("  = ", round(S.Mt.Intercept.men[i*50],5), "+ ", round(S.Mt.Slope.men[i*50],5), "* t" ))
        lines <- c(expression(paste("r"^2)), paste("            =", round(summary(current.model)$r.squared,6)) )
        text(x=c(1985,1990),y=93, labels=lines)
}

dev.off()