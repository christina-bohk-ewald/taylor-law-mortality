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

##1. Collect TL means and TL variances from observations and mortality models:  
##Note that this R source code only works after you have run "step-2-TL-Figures-1-17.R".

TL.mean.obs.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.mean.obs.men <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.obs.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.obs.men <- matrix(0,nrow=length(countrySelect),ncol=101)

TL.mean.siler.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.mean.siler.men <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.siler.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.siler.men <- matrix(0,nrow=length(countrySelect),ncol=101)

TL.mean.gompertz.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.mean.gompertz.men <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.gompertz.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.gompertz.men <- matrix(0,nrow=length(countrySelect),ncol=101)

TL.mean.gompertz.makeham.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.mean.gompertz.makeham.men <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.gompertz.makeham.wom <- matrix(0,nrow=length(countrySelect),ncol=101)
TL.var.gompertz.makeham.men <- matrix(0,nrow=length(countrySelect),ncol=101)

for(age in 1:101){ 
	for(i in 1:length(countrySelect)){
		TL.mean.obs.wom[i,age] <- TL.moments.obs.wom[[countrySelect[i]]]$mean.mx[age]
		TL.mean.obs.men[i,age] <- TL.moments.obs.men[[countrySelect[i]]]$mean.mx[age]
		TL.var.obs.wom[i,age] <- TL.moments.obs.wom[[countrySelect[i]]]$var.mx[age]
		TL.var.obs.men[i,age] <- TL.moments.obs.men[[countrySelect[i]]]$var.mx[age]
	
		TL.mean.siler.wom[i,age] <- TL.moments.siler.wom[[countrySelect[i]]]$mean.mx[age]
		TL.mean.siler.men[i,age] <- TL.moments.siler.men[[countrySelect[i]]]$mean.mx[age]
		TL.var.siler.wom[i,age] <- TL.moments.siler.wom[[countrySelect[i]]]$var.mx[age]
		TL.var.siler.men[i,age] <- TL.moments.siler.men[[countrySelect[i]]]$var.mx[age]

		TL.mean.gompertz.wom[i,age] <- TL.moments.gompertz.wom[[countrySelect[i]]]$mean.mx[age]
		TL.mean.gompertz.men[i,age] <- TL.moments.gompertz.men[[countrySelect[i]]]$mean.mx[age]
		TL.var.gompertz.wom[i,age] <- TL.moments.gompertz.wom[[countrySelect[i]]]$var.mx[age]
		TL.var.gompertz.men[i,age] <- TL.moments.gompertz.men[[countrySelect[i]]]$var.mx[age]

		TL.mean.gompertz.makeham.wom[i,age] <- TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$mean.mx[age]
		TL.mean.gompertz.makeham.men[i,age] <- TL.moments.gompertz.makeham.men[[countrySelect[i]]]$mean.mx[age]
		TL.var.gompertz.makeham.wom[i,age] <- TL.moments.gompertz.makeham.wom[[countrySelect[i]]]$var.mx[age]
		TL.var.gompertz.makeham.men[i,age] <- TL.moments.gompertz.makeham.men[[countrySelect[i]]]$var.mx[age]

	}
}

rownames(TL.mean.obs.wom) <- countrySelect
rownames(TL.var.obs.wom) <- countrySelect
rownames(TL.mean.obs.men) <- countrySelect
rownames(TL.var.obs.men) <- countrySelect

rownames(TL.mean.siler.wom) <- countrySelect
rownames(TL.var.siler.wom) <- countrySelect
rownames(TL.mean.siler.men) <- countrySelect
rownames(TL.var.siler.men) <- countrySelect

rownames(TL.mean.gompertz.wom) <- countrySelect
rownames(TL.var.gompertz.wom) <- countrySelect
rownames(TL.mean.gompertz.men) <- countrySelect
rownames(TL.var.gompertz.men) <- countrySelect

rownames(TL.mean.gompertz.makeham.wom) <- countrySelect
rownames(TL.var.gompertz.makeham.wom) <- countrySelect
rownames(TL.mean.gompertz.makeham.men) <- countrySelect
rownames(TL.var.gompertz.makeham.men) <- countrySelect

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##2. Test for Nonlinearity:

##########
#2a: Observed, women
##########

mfit.obs.men <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.obs.men) <- countrySelect
colnames(mfit.obs.men) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.obs.men[countrySelect[i],],base=10)~(log(TL.mean.obs.men[countrySelect[i],],base=10)))

	mfit.obs.men[i,1] <- coef(fit1)[1]
	mfit.obs.men[i,2] <- coef(fit1)[2]
	mfit.obs.men[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.obs.men[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.obs.men[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.obs.men[countrySelect[i],],base=10)~(log(TL.mean.obs.men[countrySelect[i],],base=10) +   I(log(TL.mean.obs.men[countrySelect[i],],base=10)^2)    ))

	mfit.obs.men[i,5] <- coef(fit2)[1]
	mfit.obs.men[i,6] <- coef(fit2)[2]
	mfit.obs.men[i,7] <- coef(fit2)[3]
	mfit.obs.men[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.obs.men[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.obs.men[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.obs.men[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.obs.men[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.obs.men[i,14] <- summary(fit2)$r.squared
	
	mfit.obs.men[i,15] <- cor(TL.var.obs.men[countrySelect[i],],TL.mean.obs.men[countrySelect[i],])
	mfit.obs.men[i,16] <- cor(log(TL.var.obs.men[countrySelect[i],],base=10),log(TL.mean.obs.men[countrySelect[i],],base=10))

	mfit.obs.men[i,17] <- cor(TL.var.obs.men[countrySelect[i],],TL.mean.obs.men[countrySelect[i],]^2)
	mfit.obs.men[i,18] <- cor(log(TL.var.obs.men[countrySelect[i],],base=10),log(TL.mean.obs.men[countrySelect[i],],base=10)^2)
}

##########
#2b: Observed, men
##########

mfit.obs.wom <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.obs.wom) <- countrySelect
colnames(mfit.obs.wom) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.obs.wom[countrySelect[i],],base=10)~(log(TL.mean.obs.wom[countrySelect[i],],base=10)))

	mfit.obs.wom[i,1] <- coef(fit1)[1]
	mfit.obs.wom[i,2] <- coef(fit1)[2]
	mfit.obs.wom[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.obs.wom[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.obs.wom[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.obs.wom[countrySelect[i],],base=10)~(log(TL.mean.obs.wom[countrySelect[i],],base=10) +   I(log(TL.mean.obs.wom[countrySelect[i],],base=10)^2)    ))

	mfit.obs.wom[i,5] <- coef(fit2)[1]
	mfit.obs.wom[i,6] <- coef(fit2)[2]
	mfit.obs.wom[i,7] <- coef(fit2)[3]
	mfit.obs.wom[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.obs.wom[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.obs.wom[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.obs.wom[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.obs.wom[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.obs.wom[i,14] <- summary(fit2)$r.squared
	
	mfit.obs.wom[i,15] <- cor(TL.var.obs.wom[countrySelect[i],],TL.mean.obs.wom[countrySelect[i],])
	mfit.obs.wom[i,16] <- cor(log(TL.var.obs.wom[countrySelect[i],],base=10),log(TL.mean.obs.wom[countrySelect[i],],base=10))

	mfit.obs.wom[i,17] <- cor(TL.var.obs.wom[countrySelect[i],],TL.mean.obs.wom[countrySelect[i],]^2)
	mfit.obs.wom[i,18] <- cor(log(TL.var.obs.wom[countrySelect[i],],base=10),log(TL.mean.obs.wom[countrySelect[i],],base=10)^2)
}

##########
#2c: Gompertz, women
##########

mfit.gompertz.men <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.gompertz.men) <- countrySelect
colnames(mfit.gompertz.men) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.gompertz.men[countrySelect[i],],base=10)~(log(TL.mean.gompertz.men[countrySelect[i],],base=10)))

	mfit.gompertz.men[i,1] <- coef(fit1)[1]
	mfit.gompertz.men[i,2] <- coef(fit1)[2]
	mfit.gompertz.men[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.gompertz.men[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.gompertz.men[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.gompertz.men[countrySelect[i],],base=10)~(log(TL.mean.gompertz.men[countrySelect[i],],base=10) +   I(log(TL.mean.gompertz.men[countrySelect[i],],base=10)^2)    ))

	mfit.gompertz.men[i,5] <- coef(fit2)[1]
	mfit.gompertz.men[i,6] <- coef(fit2)[2]
	mfit.gompertz.men[i,7] <- coef(fit2)[3]
	mfit.gompertz.men[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.gompertz.men[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.gompertz.men[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.gompertz.men[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.gompertz.men[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.gompertz.men[i,14] <- summary(fit2)$r.squared
	
	mfit.gompertz.men[i,15] <- cor(TL.var.gompertz.men[countrySelect[i],],TL.mean.gompertz.men[countrySelect[i],])
	mfit.gompertz.men[i,16] <- cor(log(TL.var.gompertz.men[countrySelect[i],],base=10),log(TL.mean.gompertz.men[countrySelect[i],],base=10))

	mfit.gompertz.men[i,17] <- cor(TL.var.gompertz.men[countrySelect[i],],TL.mean.gompertz.men[countrySelect[i],]^2)
	mfit.gompertz.men[i,18] <- cor(log(TL.var.gompertz.men[countrySelect[i],],base=10),log(TL.mean.gompertz.men[countrySelect[i],],base=10)^2)
}

##########
#2d: Gompertz, men
##########

mfit.gompertz.wom <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.gompertz.wom) <- countrySelect
colnames(mfit.gompertz.wom) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.gompertz.wom[countrySelect[i],],base=10)~(log(TL.mean.gompertz.wom[countrySelect[i],],base=10)))

	mfit.gompertz.wom[i,1] <- coef(fit1)[1]
	mfit.gompertz.wom[i,2] <- coef(fit1)[2]
	mfit.gompertz.wom[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.gompertz.wom[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.gompertz.wom[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.gompertz.wom[countrySelect[i],],base=10)~(log(TL.mean.gompertz.wom[countrySelect[i],],base=10) +   I(log(TL.mean.gompertz.wom[countrySelect[i],],base=10)^2)    ))

	mfit.gompertz.wom[i,5] <- coef(fit2)[1]
	mfit.gompertz.wom[i,6] <- coef(fit2)[2]
	mfit.gompertz.wom[i,7] <- coef(fit2)[3]
	mfit.gompertz.wom[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.gompertz.wom[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.gompertz.wom[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.gompertz.wom[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.gompertz.wom[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.gompertz.wom[i,14] <- summary(fit2)$r.squared
	
	mfit.gompertz.wom[i,15] <- cor(TL.var.gompertz.wom[countrySelect[i],],TL.mean.gompertz.wom[countrySelect[i],])
	mfit.gompertz.wom[i,16] <- cor(log(TL.var.gompertz.wom[countrySelect[i],],base=10),log(TL.mean.gompertz.wom[countrySelect[i],],base=10))

	mfit.gompertz.wom[i,17] <- cor(TL.var.gompertz.wom[countrySelect[i],],TL.mean.gompertz.wom[countrySelect[i],]^2)
	mfit.gompertz.wom[i,18] <- cor(log(TL.var.gompertz.wom[countrySelect[i],],base=10),log(TL.mean.gompertz.wom[countrySelect[i],],base=10)^2)
}

##########
#2e: Gompertz-Makeham, women
##########

mfit.gompertz.makeham.men <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.gompertz.makeham.men) <- countrySelect
colnames(mfit.gompertz.makeham.men) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.gompertz.makeham.men[countrySelect[i],],base=10)~(log(TL.mean.gompertz.makeham.men[countrySelect[i],],base=10)))

	mfit.gompertz.makeham.men[i,1] <- coef(fit1)[1]
	mfit.gompertz.makeham.men[i,2] <- coef(fit1)[2]
	mfit.gompertz.makeham.men[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.gompertz.makeham.men[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.gompertz.makeham.men[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.gompertz.makeham.men[countrySelect[i],],base=10)~(log(TL.mean.gompertz.makeham.men[countrySelect[i],],base=10) +   I(log(TL.mean.gompertz.makeham.men[countrySelect[i],],base=10)^2)    ))

	mfit.gompertz.makeham.men[i,5] <- coef(fit2)[1]
	mfit.gompertz.makeham.men[i,6] <- coef(fit2)[2]
	mfit.gompertz.makeham.men[i,7] <- coef(fit2)[3]
	mfit.gompertz.makeham.men[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.gompertz.makeham.men[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.gompertz.makeham.men[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.gompertz.makeham.men[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.gompertz.makeham.men[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.gompertz.makeham.men[i,14] <- summary(fit2)$r.squared
	
	mfit.gompertz.makeham.men[i,15] <- cor(TL.var.gompertz.makeham.men[countrySelect[i],],TL.mean.gompertz.makeham.men[countrySelect[i],])
	mfit.gompertz.makeham.men[i,16] <- cor(log(TL.var.gompertz.makeham.men[countrySelect[i],],base=10),log(TL.mean.gompertz.makeham.men[countrySelect[i],],base=10))

	mfit.gompertz.makeham.men[i,17] <- cor(TL.var.gompertz.makeham.men[countrySelect[i],],TL.mean.gompertz.makeham.men[countrySelect[i],]^2)
	mfit.gompertz.makeham.men[i,18] <- cor(log(TL.var.gompertz.makeham.men[countrySelect[i],],base=10),log(TL.mean.gompertz.makeham.men[countrySelect[i],],base=10)^2)
}

##########
#2f: Gompertz-Makeham, men
##########

mfit.gompertz.makeham.wom <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.gompertz.makeham.wom) <- countrySelect
colnames(mfit.gompertz.makeham.wom) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.gompertz.makeham.wom[countrySelect[i],],base=10)~(log(TL.mean.gompertz.makeham.wom[countrySelect[i],],base=10)))

	mfit.gompertz.makeham.wom[i,1] <- coef(fit1)[1]
	mfit.gompertz.makeham.wom[i,2] <- coef(fit1)[2]
	mfit.gompertz.makeham.wom[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.gompertz.makeham.wom[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.gompertz.makeham.wom[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.gompertz.makeham.wom[countrySelect[i],],base=10)~(log(TL.mean.gompertz.makeham.wom[countrySelect[i],],base=10) +   I(log(TL.mean.gompertz.makeham.wom[countrySelect[i],],base=10)^2)    ))

	mfit.gompertz.makeham.wom[i,5] <- coef(fit2)[1]
	mfit.gompertz.makeham.wom[i,6] <- coef(fit2)[2]
	mfit.gompertz.makeham.wom[i,7] <- coef(fit2)[3]
	mfit.gompertz.makeham.wom[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.gompertz.makeham.wom[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.gompertz.makeham.wom[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.gompertz.makeham.wom[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.gompertz.makeham.wom[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.gompertz.makeham.wom[i,14] <- summary(fit2)$r.squared
	
	mfit.gompertz.makeham.wom[i,15] <- cor(TL.var.gompertz.makeham.wom[countrySelect[i],],TL.mean.gompertz.makeham.wom[countrySelect[i],])
	mfit.gompertz.makeham.wom[i,16] <- cor(log(TL.var.gompertz.makeham.wom[countrySelect[i],],base=10),log(TL.mean.gompertz.makeham.wom[countrySelect[i],],base=10))

	mfit.gompertz.makeham.wom[i,17] <- cor(TL.var.gompertz.makeham.wom[countrySelect[i],],TL.mean.gompertz.makeham.wom[countrySelect[i],]^2)
	mfit.gompertz.makeham.wom[i,18] <- cor(log(TL.var.gompertz.makeham.wom[countrySelect[i],],base=10),log(TL.mean.gompertz.makeham.wom[countrySelect[i],],base=10)^2)
}

##########
#2g: Siler, women
##########

mfit.siler.men <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.siler.men) <- countrySelect
colnames(mfit.siler.men) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.siler.men[countrySelect[i],],base=10)~(log(TL.mean.siler.men[countrySelect[i],],base=10)))

	mfit.siler.men[i,1] <- coef(fit1)[1]
	mfit.siler.men[i,2] <- coef(fit1)[2]
	mfit.siler.men[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.siler.men[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.siler.men[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.siler.men[countrySelect[i],],base=10)~(log(TL.mean.siler.men[countrySelect[i],],base=10) +   I(log(TL.mean.siler.men[countrySelect[i],],base=10)^2)    ))

	mfit.siler.men[i,5] <- coef(fit2)[1]
	mfit.siler.men[i,6] <- coef(fit2)[2]
	mfit.siler.men[i,7] <- coef(fit2)[3]
	mfit.siler.men[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.siler.men[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.siler.men[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.siler.men[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.siler.men[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.siler.men[i,14] <- summary(fit2)$r.squared
	
	mfit.siler.men[i,15] <- cor(TL.var.siler.men[countrySelect[i],],TL.mean.siler.men[countrySelect[i],])
	mfit.siler.men[i,16] <- cor(log(TL.var.siler.men[countrySelect[i],],base=10),log(TL.mean.siler.men[countrySelect[i],],base=10))

	mfit.siler.men[i,17] <- cor(TL.var.siler.men[countrySelect[i],],TL.mean.siler.men[countrySelect[i],]^2)
	mfit.siler.men[i,18] <- cor(log(TL.var.siler.men[countrySelect[i],],base=10),log(TL.mean.siler.men[countrySelect[i],],base=10)^2)	
}

##########
#2h: Siler, men
##########

mfit.siler.wom <- matrix(0,nrow=length(countrySelect), ncol=18)
rownames(mfit.siler.wom) <- countrySelect
colnames(mfit.siler.wom) <- c("linIntercept","linBeta","pValLinInt","pValLinBeta","quadIntercept","quadBeta","quadBeta2","pValQuadInt","pValQuadBeta","pValQuadBeta2","LinAIC","QuadAIC","LinRsquared","QuadRsquared","linCorR","linlogCorR","quadCorR","quadlogCorR")

for(i in 1:length(countrySelect)){
	
	print(countrySelect[i])
	
	fit1 <- lm(log(TL.var.siler.wom[countrySelect[i],],base=10)~(log(TL.mean.siler.wom[countrySelect[i],],base=10)))

	mfit.siler.wom[i,1] <- coef(fit1)[1]
	mfit.siler.wom[i,2] <- coef(fit1)[2]
	mfit.siler.wom[i,3] <- summary(fit1)$coefficients[1,4]
	mfit.siler.wom[i,4] <- summary(fit1)$coefficients[2,4]
	mfit.siler.wom[i,13] <- summary(fit1)$r.squared
	
	fit2 <- lm(log(TL.var.siler.wom[countrySelect[i],],base=10)~(log(TL.mean.siler.wom[countrySelect[i],],base=10) +   I(log(TL.mean.siler.wom[countrySelect[i],],base=10)^2)    ))

	mfit.siler.wom[i,5] <- coef(fit2)[1]
	mfit.siler.wom[i,6] <- coef(fit2)[2]
	mfit.siler.wom[i,7] <- coef(fit2)[3]
	mfit.siler.wom[i,8] <- summary(fit2)$coefficients[1,4]
	mfit.siler.wom[i,9] <- summary(fit2)$coefficients[2,4]
	mfit.siler.wom[i,10] <- summary(fit2)$coefficients[3,4]
	
	mfit.siler.wom[i,11] <- AIC(fit1,fit2)[1,2]
	mfit.siler.wom[i,12] <- AIC(fit1,fit2)[2,2]
		
	mfit.siler.wom[i,14] <- summary(fit2)$r.squared
	
	mfit.siler.wom[i,15] <- cor(TL.var.siler.wom[countrySelect[i],],TL.mean.siler.wom[countrySelect[i],])
	mfit.siler.wom[i,16] <- cor(log(TL.var.siler.wom[countrySelect[i],],base=10),log(TL.mean.siler.wom[countrySelect[i],],base=10))

	mfit.siler.wom[i,17] <- cor(TL.var.siler.wom[countrySelect[i],],TL.mean.siler.wom[countrySelect[i],]^2)
	mfit.siler.wom[i,18] <- cor(log(TL.var.siler.wom[countrySelect[i],],base=10),log(TL.mean.siler.wom[countrySelect[i],],base=10)^2)	
}

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##3. Analysis of Covariance:

###############################
#3a: construct first dataframe: dataFrameCov:
###############################

country <- matrix(rep(c(1:12),times=1,each=101),nrow=1212,ncol=1)
age <- matrix(rep(c(0:100),times=12,each=1),nrow=1212,ncol=1)

Emx.obs.wom <- 0
Vmx.obs.wom <- 0
Emx.siler.wom <- 0
Vmx.siler.wom <- 0
Emx.gompertz.wom <- 0
Vmx.gompertz.wom <- 0
Emx.gompertz.makeham.wom <- 0
Vmx.gompertz.makeham.wom <- 0

Emx.obs.men <- 0
Vmx.obs.men <- 0
Emx.siler.men <- 0
Vmx.siler.men <- 0
Emx.gompertz.men <- 0
Vmx.gompertz.men <- 0
Emx.gompertz.makeham.men <- 0
Vmx.gompertz.makeham.men <- 0

for(i in 1:12){

	Emx.obs.wom <- c(Emx.obs.wom,t(TL.mean.obs.wom)[,i])
	Vmx.obs.wom <- c(Vmx.obs.wom,t(TL.var.obs.wom)[,i])
	Emx.siler.wom <- c(Emx.siler.wom,t(TL.mean.siler.wom)[,i])
	Vmx.siler.wom <- c(Vmx.siler.wom,t(TL.var.siler.wom)[,i])
	Emx.gompertz.wom <- c(Emx.gompertz.wom,t(TL.mean.gompertz.wom)[,i])
	Vmx.gompertz.wom <- c(Vmx.gompertz.wom,t(TL.var.gompertz.wom)[,i])
	Emx.gompertz.makeham.wom <- c(Emx.gompertz.makeham.wom,t(TL.mean.gompertz.makeham.wom)[,i])
	Vmx.gompertz.makeham.wom <- c(Vmx.gompertz.makeham.wom,t(TL.var.gompertz.makeham.wom)[,i])

	Emx.obs.men <- c(Emx.obs.men,t(TL.mean.obs.men)[,i])
	Vmx.obs.men <- c(Vmx.obs.men,t(TL.var.obs.men)[,i])
	Emx.siler.men <- c(Emx.siler.men,t(TL.mean.siler.men)[,i])
	Vmx.siler.men <- c(Vmx.siler.men,t(TL.var.siler.men)[,i])
	Emx.gompertz.men <- c(Emx.gompertz.men,t(TL.mean.gompertz.men)[,i])
	Vmx.gompertz.men <- c(Vmx.gompertz.men,t(TL.var.gompertz.men)[,i])
	Emx.gompertz.makeham.men <- c(Emx.gompertz.makeham.men,t(TL.mean.gompertz.makeham.men)[,i])
	Vmx.gompertz.makeham.men <- c(Vmx.gompertz.makeham.men,t(TL.var.gompertz.makeham.men)[,i])
}

Emx.obs.wom <- matrix(Emx.obs.wom[-1],nrow=1212,ncol=1)
Vmx.obs.wom <- matrix(Vmx.obs.wom[-1],nrow=1212,ncol=1)
Emx.siler.wom <- matrix(Emx.siler.wom[-1],nrow=1212,ncol=1)
Vmx.siler.wom <- matrix(Vmx.siler.wom[-1],nrow=1212,ncol=1)
Emx.gompertz.wom <- matrix(Emx.gompertz.wom[-1],nrow=1212,ncol=1)
Vmx.gompertz.wom <- matrix(Vmx.gompertz.wom[-1],nrow=1212,ncol=1)
Emx.gompertz.makeham.wom <- matrix(Emx.gompertz.makeham.wom[-1],nrow=1212,ncol=1)
Vmx.gompertz.makeham.wom <- matrix(Vmx.gompertz.makeham.wom[-1],nrow=1212,ncol=1)

Emx.obs.men <- matrix(Emx.obs.men[-1],nrow=1212,ncol=1)
Vmx.obs.men <- matrix(Vmx.obs.men[-1],nrow=1212,ncol=1)
Emx.siler.men <- matrix(Emx.siler.men[-1],nrow=1212,ncol=1)
Vmx.siler.men <- matrix(Vmx.siler.men[-1],nrow=1212,ncol=1)
Emx.gompertz.men <- matrix(Emx.gompertz.men[-1],nrow=1212,ncol=1)
Vmx.gompertz.men <- matrix(Vmx.gompertz.men[-1],nrow=1212,ncol=1)
Emx.gompertz.makeham.men <- matrix(Emx.gompertz.makeham.men[-1],nrow=1212,ncol=1)
Vmx.gompertz.makeham.men <- matrix(Vmx.gompertz.makeham.men[-1],nrow=1212,ncol=1)

dataCov <- cbind(country,age,Emx.obs.men,Vmx.obs.men,Emx.siler.men,Vmx.siler.men,Emx.gompertz.men,Vmx.gompertz.men,Emx.gompertz.makeham.men,Vmx.gompertz.makeham.men,
				Emx.obs.wom,Vmx.obs.wom,Emx.siler.wom,Vmx.siler.wom,Emx.gompertz.wom,Vmx.gompertz.wom,Emx.gompertz.makeham.wom,Vmx.gompertz.makeham.wom)
colnames(dataCov) <- c("country","age","Emx.obs.men","Vmx.obs.men","Emx.siler.men","Vmx.siler.men","Emx.gompertz.men","Vmx.gompertz.men","Emx.gompertz.makeham.men","Vmx.gompertz.makeham.men",
						"Emx.obs.wom","Vmx.obs.wom","Emx.siler.wom","Vmx.siler.wom","Emx.gompertz.wom","Vmx.gompertz.wom","Emx.gompertz.makeham.wom","Vmx.gompertz.makeham.wom")

dataFrameCov <- as.data.frame(dataCov)
attach(dataFrameCov)

##set country as factor in lm:

dataFrameCov$country = factor(as.numeric(dataFrameCov$country))
dataFrameCov$age = factor(as.numeric(dataFrameCov$age))

##reduce dataFrameCov to selected countries or ages:

agesSelect <- c(1, (1:10)*10+1)

dfCov.SelectAges <- dataFrameCov[dataFrameCov$age %in% agesSelect,]
dfCov.SelectAges$country = factor(as.numeric(dfCov.SelectAges$country))
dfCov.SelectAges$age = factor(as.numeric(dfCov.SelectAges$age))

dfCov.SelectCountries <- dataFrameCov[dataFrameCov$country %in% match(countrySelect,countrySelect),]
dfCov.SelectCountries$country = factor(as.numeric(dfCov.SelectCountries$country) ) 
dfCov.SelectCountries$age = factor(as.numeric(dfCov.SelectCountries$age))

###############################
#3b: construct (final) dataframe: dataFrameCov2:
###############################

##Differences between models and / or sexes?

country <- matrix(rep(c(1:12),times=1,each=101),nrow=1212,ncol=1)
age <- matrix(rep(c(0:100),times=12,each=1),nrow=1212,ncol=1)
sex1 <- matrix(rep(c(1),times=1212),nrow=1212,ncol=1)
sex2 <- matrix(rep(c(2),times=1212),nrow=1212,ncol=1)
model0 <- matrix(rep(c(0),times=1212),nrow=1212,ncol=1) 
model1 <- matrix(rep(c(1),times=1212),nrow=1212,ncol=1) 
model2 <- matrix(rep(c(2),times=1212),nrow=1212,ncol=1) 
model3 <- matrix(rep(c(3),times=1212),nrow=1212,ncol=1) 

df.Emx.obs.wom <- cbind(country,age,sex1,model0,Emx.obs.wom,Vmx.obs.wom)
df.Emx.obs.men <- cbind(country,age,sex2,model0,Emx.obs.men,Vmx.obs.men)
df.Emx.siler.wom <- cbind(country,age,sex1,model1,Emx.siler.wom,Vmx.siler.wom)
df.Emx.siler.men <- cbind(country,age,sex2,model1,Emx.siler.men,Vmx.siler.men)
df.Emx.gompertz.wom <- cbind(country,age,sex1,model2,Emx.gompertz.wom,Vmx.gompertz.wom)
df.Emx.gompertz.men <- cbind(country,age,sex2,model2,Emx.gompertz.men,Vmx.gompertz.men)
df.Emx.gompertz.makeham.wom <- cbind(country,age,sex1,model3,Emx.gompertz.makeham.wom,Vmx.gompertz.makeham.wom)
df.Emx.gompertz.makeham.men <- cbind(country,age,sex2,model3,Emx.gompertz.makeham.men,Vmx.gompertz.makeham.men)

dataCov2 <- rbind(df.Emx.obs.wom,df.Emx.obs.men,
					df.Emx.siler.wom,df.Emx.siler.men,
					df.Emx.gompertz.wom,df.Emx.gompertz.men,
					df.Emx.gompertz.makeham.wom,df.Emx.gompertz.makeham.men)
colnames(dataCov2) <- c("country","age","sex","model","Emx","Vmx")

dataFrameCov2 <- as.data.frame(dataCov2)
attach(dataFrameCov2)

##set country as factor in lm:

dataFrameCov2$country = factor(as.numeric(dataFrameCov2$country))
dataFrameCov2$age = factor(as.numeric(dataFrameCov2$age))
dataFrameCov2$sex = factor(as.numeric(dataFrameCov2$sex))
dataFrameCov2$model = factor(as.numeric(dataFrameCov2$model))

#################################
#3c: Regression with interaction:
#################################

###############
###############
##Table 1 of Cohen, Bohk-Ewald, Rau (2018):
###############
###############

##Do the slopes of TL differ between observations and fitted mortality models?

####
####
####Table 1, Women:
####
####

#Model differences by sex (=Women) and country (=all countries) (Observed as reference):

fit.model.country.all.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.all.wom)

#Model differences by sex (=Women) and country (=DNK) (Observed as reference):

fit.model.country.dnk1.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==1&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==1&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.dnk1.wom)

#Model differences by sex (=Women) and country (=FRA) (Observed as reference):

fit.model.country.fra2.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==2&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==2&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.fra2.wom) 

#Model differences by sex (=Women) and country (=GDR) (Observed as reference):

fit.model.country.gdr3.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==3&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==3&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.gdr3.wom)

#Model differences by sex (=Women) and country (=FRG) (Observed as reference):

fit.model.country.frg4.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==4&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==4&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.frg4.wom)

#Model differences by sex (=Women) and country (=HUN) (Observed as reference):

fit.model.country.hun5.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==5&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==5&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.hun5.wom)

#Model differences by sex (=Women) and country (=ITA) (Observed as reference):

fit.model.country.ita6.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==6&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==6&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.ita6.wom)

#Model differences by sex (=Women) and country (=JPN) (Observed as reference):

fit.model.country.jpn7.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==7&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==7&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.jpn7.wom)

#Model differences by sex (=Women) and country (=POL) (Observed as reference):

fit.model.country.pol8.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==8&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==8&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.pol8.wom)

#Model differences by sex (=Women) and country (=RUS) (Observed as reference):

fit.model.country.rus9.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==9&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==9&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.rus9.wom)

#Model differences by sex (=Women) and country (=SWE) (Observed as reference):

fit.model.country.swe10.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==10&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==10&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.swe10.wom)

#Model differences by sex (=Women) and country (=UK) (Observed as reference):

fit.model.country.uk11.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==11&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==11&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.uk11.wom)

#Model differences by sex (=Women) and country (=USA) (Observed as reference):

fit.model.country.usa12.wom <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$sex==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$sex==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==12&dataFrameCov2$sex==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$sex==1],base=10):dataFrameCov2$model[dataFrameCov2$country==12&dataFrameCov2$sex==1], 
data=dataFrameCov2)

summary(fit.model.country.usa12.wom)

####
####
####Table 1, Men:
####
####

#Model differences by sex (=Men ) and country (=all countries) (Observed as reference):

fit.model.country.all.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.all.men)

#Model differences by sex (=Men ) and country (=DNK) (Observed as reference):

fit.model.country.dnk1.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==1&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==1&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.dnk1.men)

#Model differences by sex (=Men ) and country (=FRA) (Observed as reference):

fit.model.country.fra2.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==2&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==2&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.fra2.men) 

#Model differences by sex (=Men ) and country (=GDR) (Observed as reference):

fit.model.country.gdr3.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==3&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==3&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.gdr3.men)

#Model differences by sex (=Men ) and country (=FRG) (Observed as reference):

fit.model.country.frg4.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==4&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==4&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.frg4.men)

#Model differences by sex (=Men ) and country (=HUN) (Observed as reference):

fit.model.country.hun5.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==5&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==5&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.hun5.men)

#Model differences by sex (=Men ) and country (=ITA) (Observed as reference):

fit.model.country.ita6.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==6&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==6&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.ita6.men)

#Model differences by sex (=Men ) and country (=JPN) (Observed as reference):

fit.model.country.jpn7.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==7&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==7&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.jpn7.men)

#Model differences by sex (=Men ) and country (=POL) (Observed as reference):

fit.model.country.pol8.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==8&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==8&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.pol8.men)

#Model differences by sex (=Men ) and country (=RUS) (Observed as reference):

fit.model.country.rus9.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==9&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==9&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.rus9.men)

#Model differences by sex (=Men ) and country (=SWE) (Observed as reference):

fit.model.country.swe10.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==10&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==10&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.swe10.men)

#Model differences by sex (=Men ) and country (=UK) (Observed as reference):

fit.model.country.uk11.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==11&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==11&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.uk11.men)

#Model differences by sex (=Men ) and country (=USA) (Observed as reference):

fit.model.country.usa12.men <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$sex==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$sex==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==12&dataFrameCov2$sex==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$sex==2],base=10):dataFrameCov2$model[dataFrameCov2$country==12&dataFrameCov2$sex==2], 
data=dataFrameCov2)

summary(fit.model.country.usa12.men)

####
####
####Table 1, Both sexes:
####
####

#Model differences by sex (=Both ) and country (=all countries) (Observed as reference):

fit.model.country.all.both <- lm(log(dataFrameCov2$Vmx,base=10) ~ log(dataFrameCov2$Emx,base=10) + dataFrameCov2$model + log(dataFrameCov2$Emx,base=10):dataFrameCov2$model, 
data=dataFrameCov2)

summary(fit.model.country.all.both)

#Model differences by sex (=Both ) and country (=DNK) (Observed as reference):

fit.model.country.dnk1.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10) + dataFrameCov2$model[dataFrameCov2$country==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10):dataFrameCov2$model[dataFrameCov2$country==1], 
data=dataFrameCov2)

summary(fit.model.country.dnk1.both)

#Model differences by sex (=Both ) and country (=FRA) (Observed as reference):

fit.model.country.fra2.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10) + dataFrameCov2$model[dataFrameCov2$country==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10):dataFrameCov2$model[dataFrameCov2$country==2], 
data=dataFrameCov2)

summary(fit.model.country.fra2.both) 

#Model differences by sex (=Both ) and country (=GDR) (Observed as reference):

fit.model.country.gdr3.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10) + dataFrameCov2$model[dataFrameCov2$country==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10):dataFrameCov2$model[dataFrameCov2$country==3], 
data=dataFrameCov2)

summary(fit.model.country.gdr3.both)

#Model differences by sex (=Both ) and country (=FRG) (Observed as reference):

fit.model.country.frg4.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10) + dataFrameCov2$model[dataFrameCov2$country==4] + log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10):dataFrameCov2$model[dataFrameCov2$country==4], 
data=dataFrameCov2)

summary(fit.model.country.frg4.both)

#Model differences by sex (=Both ) and country (=HUN) (Observed as reference):

fit.model.country.hun5.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10) + dataFrameCov2$model[dataFrameCov2$country==5] + log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10):dataFrameCov2$model[dataFrameCov2$country==5], 
data=dataFrameCov2)

summary(fit.model.country.hun5.both)

#Model differences by sex (=Both ) and country (=ITA) (Observed as reference):

fit.model.country.ita6.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10) + dataFrameCov2$model[dataFrameCov2$country==6] + log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10):dataFrameCov2$model[dataFrameCov2$country==6], 
data=dataFrameCov2)

summary(fit.model.country.ita6.both)

#Model differences by sex (=Both ) and country (=JPN) (Observed as reference):

fit.model.country.jpn7.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10) + dataFrameCov2$model[dataFrameCov2$country==7] + log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10):dataFrameCov2$model[dataFrameCov2$country==7], 
data=dataFrameCov2)

summary(fit.model.country.jpn7.both)

#Model differences by sex (=Both ) and country (=POL) (Observed as reference):

fit.model.country.pol8.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10) + dataFrameCov2$model[dataFrameCov2$country==8] + log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10):dataFrameCov2$model[dataFrameCov2$country==8], 
data=dataFrameCov2)

summary(fit.model.country.pol8.both)

#Model differences by sex (=Both ) and country (=RUS) (Observed as reference):

fit.model.country.rus9.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10) + dataFrameCov2$model[dataFrameCov2$country==9] + log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10):dataFrameCov2$model[dataFrameCov2$country==9], 
data=dataFrameCov2)

summary(fit.model.country.rus9.both)

#Model differences by sex (=Both ) and country (=SWE) (Observed as reference):

fit.model.country.swe10.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10) + dataFrameCov2$model[dataFrameCov2$country==10] + log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10):dataFrameCov2$model[dataFrameCov2$country==10], 
data=dataFrameCov2)

summary(fit.model.country.swe10.both)

#Model differences by sex (=Both ) and country (=UK) (Observed as reference):

fit.model.country.uk11.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10) + dataFrameCov2$model[dataFrameCov2$country==11] + log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10):dataFrameCov2$model[dataFrameCov2$country==11], 
data=dataFrameCov2)

summary(fit.model.country.uk11.both)

#Model differences by sex (=Both ) and country (=USA) (Observed as reference):

fit.model.country.usa12.both <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10) + dataFrameCov2$model[dataFrameCov2$country==12] + log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10):dataFrameCov2$model[dataFrameCov2$country==12], 
data=dataFrameCov2)

summary(fit.model.country.usa12.both)

################
################
##Table 2 of Cohen, Bohk-Ewald, Rau (2018):
################
################

##Do the slopes of TL differ between women and men?

####
####
####Table 2, Observed data:
####
####

#Differences in sex by model (=observed) (Women as reference)---regardless of country:

fit.sex.model.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.model.obs)  

#Differences in sex by country (=DNK) (Women as reference)---regard model = observed:

fit.sex.country.dnk1.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.dnk1.obs)

#Differences in sex by country (=FRA) (Women as reference)---regard model = observed:

fit.sex.country.fra2.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.fra2.obs) 

#Differences in sex by country (=GDR) (Women as reference)---regard model = observed:

fit.sex.country.gdr3.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.gdr3.obs)

#Differences in sex by country (=FRG) (Women as reference)---regard model = observed:

fit.sex.country.frg4.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.frg4.obs)

#Differences in sex by country (=HUN) (Women as reference)---regard model = observed:

fit.sex.country.hun5.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.hun5.obs)

#Differences in sex by country (=ITA) (Women as reference)---regard model = observed:

fit.sex.country.ita6.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.ita6.obs)

#Differences in sex by country (=JPN) (Women as reference)---regard model = observed:

fit.sex.country.jpn7.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.jpn7.obs)

#Differences in sex by country (=POL) (Women as reference)---regard model = observed:

fit.sex.country.pol8.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.pol8.obs)

#Differences in sex by country (=RUS) (Women as reference)---regard model = observed:

fit.sex.country.rus9.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.rus9.obs)

#Differences in sex by country (=SWE) (Women as reference)---regard model = observed:

fit.sex.country.swe10.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.swe10.obs)

#Differences in sex by country (=UK) (Women as reference)---regard model = observed:

fit.sex.country.uk11.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.uk11.obs)

#Differences in sex by country (=USA) (Women as reference)---regard model = observed:

fit.sex.country.usa12.obs <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$model==0],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==0],base=10) + dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==0] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==0],base=10):dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==0], 
data=dataFrameCov2)

summary(fit.sex.country.usa12.obs)

####
####
####Table 2, Gompertz model:
####
####

#Differences in sex by model (=gompertz) (Women as reference)---regardless of country:

fit.sex.model.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.model.gompertz)   ###differences in TL in gompertz model between women and men!!!! 

#Differences in sex by country (=DNK) (Women as reference)---regard model = gompertz:

fit.sex.country.dnk1.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.dnk1.gompertz)

#Differences in sex by country (=FRA) (Women as reference)---regard model = gompertz:

fit.sex.country.fra2.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.fra2.gompertz) 

#Differences in sex by country (=GDR) (Women as reference)---regard model = gompertz:

fit.sex.country.gdr3.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.gdr3.gompertz)

#Differences in sex by country (=FRG) (Women as reference)---regard model = gompertz:

fit.sex.country.frg4.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.frg4.gompertz)

#Differences in sex by country (=HUN) (Women as reference)---regard model = gompertz:

fit.sex.country.hun5.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.hun5.gompertz)

#Differences in sex by country (=ITA) (Women as reference)---regard model = gompertz:

fit.sex.country.ita6.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.ita6.gompertz)

#Differences in sex by country (=JPN) (Women as reference)---regard model = gompertz:

fit.sex.country.jpn7.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.jpn7.gompertz)

#Differences in sex by country (=POL) (Women as reference)---regard model = gompertz:

fit.sex.country.pol8.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.pol8.gompertz)

#Differences in sex by country (=RUS) (Women as reference)---regard model = gompertz:

fit.sex.country.rus9.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.rus9.gompertz)

#Differences in sex by country (=SWE) (Women as reference)---regard model = gompertz:

fit.sex.country.swe10.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.swe10.gompertz)

#Differences in sex by country (=UK) (Women as reference)---regard model = gompertz:

fit.sex.country.uk11.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.uk11.gompertz)

#Differences in sex by country (=USA) (Women as reference)---regard model = gompertz:

fit.sex.country.usa12.gompertz <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$model==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==2] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==2], 
data=dataFrameCov2)

summary(fit.sex.country.usa12.gompertz)

####
####
####Table 2, Gompertz-Makeham model:
####
####

#Differences in sex by model (=gompertz.makeham) (Women as reference)---regardless of country:

fit.sex.model.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.model.gompertz.makeham)

#Differences in sex by country (=DNK) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.dnk1.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.dnk1.gompertz.makeham)

#Differences in sex by country (=FRA) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.fra2.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.fra2.gompertz.makeham) 

#Differences in sex by country (=GDR) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.gdr3.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.gdr3.gompertz.makeham)

#Differences in sex by country (=FRG) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.frg4.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.frg4.gompertz.makeham)

#Differences in sex by country (=HUN) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.hun5.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.hun5.gompertz.makeham)

#Differences in sex by country (=ITA) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.ita6.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.ita6.gompertz.makeham)

#Differences in sex by country (=JPN) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.jpn7.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.jpn7.gompertz.makeham)

#Differences in sex by country (=POL) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.pol8.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.pol8.gompertz.makeham)

#Differences in sex by country (=RUS) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.rus9.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.rus9.gompertz.makeham)

#Differences in sex by country (=SWE) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.swe10.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.swe10.gompertz.makeham)

#Differences in sex by country (=UK) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.uk11.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.uk11.gompertz.makeham)

#Differences in sex by country (=USA) (Women as reference)---regard model = gompertz.makeham:

fit.sex.country.usa12.gompertz.makeham <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$model==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==3] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==3], 
data=dataFrameCov2)

summary(fit.sex.country.usa12.gompertz.makeham)

####
####
####Table 2, Siler model:
####
####

#Differences in sex by model (=siler) (Women as reference)---regardless of country:

fit.sex.model.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.model.siler)   

#Differences in sex by country (=DNK) (Women as reference)---regard model = siler:

fit.sex.country.dnk1.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==1&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==1&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.dnk1.siler)

#Differences in sex by country (=FRA) (Women as reference)---regard model = siler:

fit.sex.country.fra2.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==2&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==2&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.fra2.siler) 

#Differences in sex by country (=GDR) (Women as reference)---regard model = siler:

fit.sex.country.gdr3.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==3&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==3&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.gdr3.siler)

#Differences in sex by country (=FRG) (Women as reference)---regard model = siler:

fit.sex.country.frg4.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==4&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==4&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.frg4.siler)

#Differences in sex by country (=HUN) (Women as reference)---regard model = siler:

fit.sex.country.hun5.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==5&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==5&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.hun5.siler)

#Differences in sex by country (=ITA) (Women as reference)---regard model = siler:

fit.sex.country.ita6.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==6&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==6&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.ita6.siler)

#Differences in sex by country (=JPN) (Women as reference)---regard model = siler:

fit.sex.country.jpn7.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==7&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==7&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.jpn7.siler)

#Differences in sex by country (=POL) (Women as reference)---regard model = siler:

fit.sex.country.pol8.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==8&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==8&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.pol8.siler)

#Differences in sex by country (=RUS) (Women as reference)---regard model = siler:

fit.sex.country.rus9.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==9&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==9&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.rus9.siler)

#Differences in sex by country (=SWE) (Women as reference)---regard model = siler:

fit.sex.country.swe10.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==10&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==10&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.swe10.siler)

#Differences in sex by country (=UK) (Women as reference)---regard model = siler:

fit.sex.country.uk11.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==11&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==11&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.uk11.siler)

#Differences in sex by country (=USA) (Women as reference)---regard model = siler:

fit.sex.country.usa12.siler <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12&dataFrameCov2$model==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==1] + log(dataFrameCov2$Emx[dataFrameCov2$country==12&dataFrameCov2$model==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==12&dataFrameCov2$model==1], 
data=dataFrameCov2)

summary(fit.sex.country.usa12.siler)

################
################
##Table 3 of Cohen, Bohk-Ewald, Rau (2018):
################
################

##Do the sex differences in slopes of TL differ between observations and the mortality models?

####
####
####Table 3, Siler (=model 1), Gompertz (=model 2), Gompertz-Makeham (=model 3):
####
####


#Differences in sex by country (=all countries) (Women as reference)---reference model = observed:


fit.sex.country.all.obs.SM <- lm(log(dataFrameCov2$Vmx,base=10) ~ log(dataFrameCov2$Emx,base=10) + dataFrameCov2$sex + dataFrameCov2$model + 
									log(dataFrameCov2$Emx,base=10):dataFrameCov2$sex +
									log(dataFrameCov2$Emx,base=10):dataFrameCov2$model +
	 								dataFrameCov2$sex:dataFrameCov2$model +
	 								log(dataFrameCov2$Emx,base=10):dataFrameCov2$sex:dataFrameCov2$model, 
									data=dataFrameCov2)

summary(fit.sex.country.all.obs.SM)

#Differences in sex by country (=DNK) (Women as reference)---reference model = observed:

fit.sex.country.dnk1.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==1],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10) + dataFrameCov2$sex[dataFrameCov2$country==1] + dataFrameCov2$model[dataFrameCov2$country==1] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==1] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10):dataFrameCov2$model[dataFrameCov2$country==1] +
	 								dataFrameCov2$sex[dataFrameCov2$country==1]:dataFrameCov2$model[dataFrameCov2$country==1] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==1],base=10):dataFrameCov2$sex[dataFrameCov2$country==1]:dataFrameCov2$model[dataFrameCov2$country==1], 
									data=dataFrameCov2)

summary(fit.sex.country.dnk1.obs.SM)

#Differences in sex by country (=FRA) (Women as reference)---reference model = observed:

fit.sex.country.fra2.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==2],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10) + dataFrameCov2$sex[dataFrameCov2$country==2] + dataFrameCov2$model[dataFrameCov2$country==2] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==2] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10):dataFrameCov2$model[dataFrameCov2$country==2] +
	 								dataFrameCov2$sex[dataFrameCov2$country==2]:dataFrameCov2$model[dataFrameCov2$country==2] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==2],base=10):dataFrameCov2$sex[dataFrameCov2$country==2]:dataFrameCov2$model[dataFrameCov2$country==2], 
									data=dataFrameCov2)

summary(fit.sex.country.fra2.obs.SM)

#Differences in sex by country (=GDR) (Women as reference)---reference model = observed:

fit.sex.country.gdr3.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==3],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10) + dataFrameCov2$sex[dataFrameCov2$country==3] + dataFrameCov2$model[dataFrameCov2$country==3] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==3] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10):dataFrameCov2$model[dataFrameCov2$country==3] +
	 								dataFrameCov2$sex[dataFrameCov2$country==3]:dataFrameCov2$model[dataFrameCov2$country==3] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==3],base=10):dataFrameCov2$sex[dataFrameCov2$country==3]:dataFrameCov2$model[dataFrameCov2$country==3], 
									data=dataFrameCov2)

summary(fit.sex.country.gdr3.obs.SM)

#Differences in sex by country (=FRG) (Women as reference)---reference model = observed:

fit.sex.country.frg4.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==4],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10) + dataFrameCov2$sex[dataFrameCov2$country==4] + dataFrameCov2$model[dataFrameCov2$country==4] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10):dataFrameCov2$sex[dataFrameCov2$country==4] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10):dataFrameCov2$model[dataFrameCov2$country==4] +
	 								dataFrameCov2$sex[dataFrameCov2$country==4]:dataFrameCov2$model[dataFrameCov2$country==4] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==4],base=10):dataFrameCov2$sex[dataFrameCov2$country==4]:dataFrameCov2$model[dataFrameCov2$country==4], 
									data=dataFrameCov2)

summary(fit.sex.country.frg4.obs.SM)

#Differences in sex by country (=HUN) (Women as reference)---reference model = observed:

fit.sex.country.hun5.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==5],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10) + dataFrameCov2$sex[dataFrameCov2$country==5] + dataFrameCov2$model[dataFrameCov2$country==5] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10):dataFrameCov2$sex[dataFrameCov2$country==5] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10):dataFrameCov2$model[dataFrameCov2$country==5] +
	 								dataFrameCov2$sex[dataFrameCov2$country==5]:dataFrameCov2$model[dataFrameCov2$country==5] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==5],base=10):dataFrameCov2$sex[dataFrameCov2$country==5]:dataFrameCov2$model[dataFrameCov2$country==5], 
									data=dataFrameCov2)

summary(fit.sex.country.hun5.obs.SM)

#Differences in sex by country (=ITA) (Women as reference)---reference model = observed:

fit.sex.country.ita6.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==6],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10) + dataFrameCov2$sex[dataFrameCov2$country==6] + dataFrameCov2$model[dataFrameCov2$country==6] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10):dataFrameCov2$sex[dataFrameCov2$country==6] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10):dataFrameCov2$model[dataFrameCov2$country==6] +
	 								dataFrameCov2$sex[dataFrameCov2$country==6]:dataFrameCov2$model[dataFrameCov2$country==6] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==6],base=10):dataFrameCov2$sex[dataFrameCov2$country==6]:dataFrameCov2$model[dataFrameCov2$country==6], 
									data=dataFrameCov2)

summary(fit.sex.country.ita6.obs.SM)

#Differences in sex by country (=JPN) (Women as reference)---reference model = observed:

fit.sex.country.jpn7.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==7],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10) + dataFrameCov2$sex[dataFrameCov2$country==7] + dataFrameCov2$model[dataFrameCov2$country==7] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10):dataFrameCov2$sex[dataFrameCov2$country==7] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10):dataFrameCov2$model[dataFrameCov2$country==7] +
	 								dataFrameCov2$sex[dataFrameCov2$country==7]:dataFrameCov2$model[dataFrameCov2$country==7] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==7],base=10):dataFrameCov2$sex[dataFrameCov2$country==7]:dataFrameCov2$model[dataFrameCov2$country==7], 
									data=dataFrameCov2)

summary(fit.sex.country.jpn7.obs.SM)

#Differences in sex by country (=POL) (Women as reference)---reference model = observed:

fit.sex.country.pol8.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==8],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10) + dataFrameCov2$sex[dataFrameCov2$country==8] + dataFrameCov2$model[dataFrameCov2$country==8] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10):dataFrameCov2$sex[dataFrameCov2$country==8] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10):dataFrameCov2$model[dataFrameCov2$country==8] +
	 								dataFrameCov2$sex[dataFrameCov2$country==8]:dataFrameCov2$model[dataFrameCov2$country==8] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==8],base=10):dataFrameCov2$sex[dataFrameCov2$country==8]:dataFrameCov2$model[dataFrameCov2$country==8], 
									data=dataFrameCov2)

summary(fit.sex.country.pol8.obs.SM)

#Differences in sex by country (=RUS) (Women as reference)---reference model = observed:

fit.sex.country.rus9.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==9],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10) + dataFrameCov2$sex[dataFrameCov2$country==9] + dataFrameCov2$model[dataFrameCov2$country==9] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10):dataFrameCov2$sex[dataFrameCov2$country==9] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10):dataFrameCov2$model[dataFrameCov2$country==9] +
	 								dataFrameCov2$sex[dataFrameCov2$country==9]:dataFrameCov2$model[dataFrameCov2$country==9] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==9],base=10):dataFrameCov2$sex[dataFrameCov2$country==9]:dataFrameCov2$model[dataFrameCov2$country==9], 
									data=dataFrameCov2)

summary(fit.sex.country.rus9.obs.SM)

#Differences in sex by country (=SWE) (Women as reference)---reference model = observed:

fit.sex.country.swe10.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==10],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10) + dataFrameCov2$sex[dataFrameCov2$country==10] + dataFrameCov2$model[dataFrameCov2$country==10] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10):dataFrameCov2$sex[dataFrameCov2$country==10] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10):dataFrameCov2$model[dataFrameCov2$country==10] +
	 								dataFrameCov2$sex[dataFrameCov2$country==10]:dataFrameCov2$model[dataFrameCov2$country==10] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==10],base=10):dataFrameCov2$sex[dataFrameCov2$country==10]:dataFrameCov2$model[dataFrameCov2$country==10], 
									data=dataFrameCov2)

summary(fit.sex.country.swe10.obs.SM)

#Differences in sex by country (=UK) (Women as reference)---reference model = observed:

fit.sex.country.uk11.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==11],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10) + dataFrameCov2$sex[dataFrameCov2$country==11] + dataFrameCov2$model[dataFrameCov2$country==11] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10):dataFrameCov2$sex[dataFrameCov2$country==11] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10):dataFrameCov2$model[dataFrameCov2$country==11] +
	 								dataFrameCov2$sex[dataFrameCov2$country==11]:dataFrameCov2$model[dataFrameCov2$country==11] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==11],base=10):dataFrameCov2$sex[dataFrameCov2$country==11]:dataFrameCov2$model[dataFrameCov2$country==11], 
									data=dataFrameCov2)

summary(fit.sex.country.uk11.obs.SM)

#Differences in sex by country (=USA) (Women as reference)---reference model = observed:

fit.sex.country.usa12.obs.SM <- lm(log(dataFrameCov2$Vmx[dataFrameCov2$country==12],base=10) ~ log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10) + dataFrameCov2$sex[dataFrameCov2$country==12] + dataFrameCov2$model[dataFrameCov2$country==12] + 
									log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10):dataFrameCov2$sex[dataFrameCov2$country==12] +
									log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10):dataFrameCov2$model[dataFrameCov2$country==12] +
	 								dataFrameCov2$sex[dataFrameCov2$country==12]:dataFrameCov2$model[dataFrameCov2$country==12] +
	 								log(dataFrameCov2$Emx[dataFrameCov2$country==12],base=10):dataFrameCov2$sex[dataFrameCov2$country==12]:dataFrameCov2$model[dataFrameCov2$country==12], 
									data=dataFrameCov2)

summary(fit.sex.country.usa12.obs.SM)

#######################################################################
#######################################################################

