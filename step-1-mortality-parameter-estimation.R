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

##1. Load required data:

##Please note: to make use of the R code used in Cohen, Bohk-Ewald, Rau (2018) you will need to prepare mortality data yourself. 
##You will need four list objects: deaths.wom, deaths.men, exposures.wom, and exposures.men.
##Each list object contains the deaths (or exposures) by single years of age (here: 0 through 110) and single calendar years (here: 1950 through 2015) for women (or men) for a list of selected countries.  
##In Cohen, Bohk-Ewald, Rau (2018), we have used data of the Human Mortality Database that is available at: www.mortality.org.

deaths.wom <- list()
deaths.men <- list()
exposures.wom <- list()
exposures.men <- list()

countrySelect <- c("DNK","FRA","GDR","FRG","HUN","ITA","JPN","POL","RUS","SWE","GBR","USA")
countrySelect2 <- c("Denmark","France","East Germany","West Germany","Hungary","Italy","Japan","Poland","Russia","Sweden","United Kingdom","USA")

##
##PUT IN HERE YOUR R CODE TO PREPARE THE FOUR LIST OBJECTS: deaths.wom, deaths.men, exposures.wom, and exposures.men, AND STORE THEM IN ./dumps!
##

setwd(the.script.path)
setwd(the.dump.path)

dump(list="deaths.wom",file="deaths_wom.R")
dump(list="deaths.men",file="deaths_men.R")
dump(list="exposures.wom",file="exposures_wom.R")
dump(list="exposures.men",file="exposures_men.R")

source("deaths_wom.R")
source("deaths_men.R")
source("exposures_wom.R")
source("exposures_men.R")

#######################################################################
#######################################################################


#######################################################################
#######################################################################

##2. Load required functions and estimate parameters of the mortality models used in Cohen, Bohk-Ewald, Rau (2018):

require(DEoptim)
#install.packages("DEoptim")

setwd(the.script.path)
source("fctn-gompertz.R")
source("fctn-gompertz-makeham.R")
source("fctn-siler.R")

###########
##2a: Gompertz, women:
###########

gompertz.par.wom <- list()
gompertz.est.wom <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.Dx <- deaths.wom[[countrySelect[i]]]
	current.Nx <- exposures.wom[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=2)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.gompertz,lower=c(50,0.01),
						upper=c(110,4), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("M","b")
	rownames(current.par) <- colnames(current.Dx)

	gompertz.par.wom[[countrySelect[i]]] <- current.par


	current.gompertz.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.gompertz.estimates[,year] <- gompertz(M=current.par[year,1],
											b=current.par[year,2],
											x=current.age)
	}

	rownames(current.gompertz.estimates) <- current.age
	colnames(current.gompertz.estimates) <- rownames(current.par)

	gompertz.est.wom[[countrySelect[i]]] <- current.gompertz.estimates

}

setwd(the.dump.path)

dump(list="gompertz.par.wom",file="gompertz_par_wom.R")
dump(list="gompertz.est.wom",file="gompertz_est_wom.R")

###########
##2b: Gompertz, men:
###########

gompertz.par.men <- list()
gompertz.est.men <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.Dx <- deaths.men[[countrySelect[i]]]
	current.Nx <- exposures.men[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=2)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.gompertz,lower=c(50,0.01),
						upper=c(110,4), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("M","b")
	rownames(current.par) <- colnames(current.Dx)

	gompertz.par.men[[countrySelect[i]]] <- current.par


	current.gompertz.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.gompertz.estimates[,year] <- gompertz(M=current.par[year,1],
											b=current.par[year,2],
											x=current.age)
	}

	rownames(current.gompertz.estimates) <- current.age
	colnames(current.gompertz.estimates) <- rownames(current.par)

	gompertz.est.men[[countrySelect[i]]] <- current.gompertz.estimates

}

setwd(the.dump.path)

dump(list="gompertz.par.men",file="gompertz_par_men.R")
dump(list="gompertz.est.men",file="gompertz_est_men.R")

###########
##2c: Gompertz-Makeham, women:
###########

gompertz.makeham.par.wom <- list()
gompertz.makeham.est.wom <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.Dx <- deaths.wom[[countrySelect[i]]]
	current.Nx <- exposures.wom[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=3)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.gompertz.makeham,lower=c(0,50,0.04),
						upper=c(0.1,110,4), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("c","M","b")
	rownames(current.par) <- colnames(current.Dx)

	gompertz.makeham.par.wom[[countrySelect[i]]] <- current.par


	current.gompertz.makeham.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.gompertz.makeham.estimates[,year] <- gompertz.makeham(c=current.par[year,1],
											M=current.par[year,2],
											b=current.par[year,3],
											x=current.age)
	}

	rownames(current.gompertz.makeham.estimates) <- current.age
	colnames(current.gompertz.makeham.estimates) <- rownames(current.par)

	gompertz.makeham.est.wom[[countrySelect[i]]] <- current.gompertz.makeham.estimates

}

setwd(the.dump.path)

dump(list="gompertz.makeham.par.wom",file="gompertz_makeham_par_wom.R")
dump(list="gompertz.makeham.est.wom",file="gompertz_makeham_est_wom.R")

###########
##2d: Gompertz-Makeham, men:
###########

gompertz.makeham.par.men <- list()
gompertz.makeham.est.men <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.Dx <- deaths.men[[countrySelect[i]]]
	current.Nx <- exposures.men[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=3)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.gompertz.makeham,lower=c(0,50,0.04),
						upper=c(0.1,110,4), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("c","M","b")
	rownames(current.par) <- colnames(current.Dx)

	gompertz.makeham.par.men[[countrySelect[i]]] <- current.par


	current.gompertz.makeham.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.gompertz.makeham.estimates[,year] <- gompertz.makeham(c=current.par[year,1],
											M=current.par[year,2],
											b=current.par[year,3],
											x=current.age)
	}

	rownames(current.gompertz.makeham.estimates) <- current.age
	colnames(current.gompertz.makeham.estimates) <- rownames(current.par)

	gompertz.makeham.est.men[[countrySelect[i]]] <- current.gompertz.makeham.estimates

}

setwd(the.dump.path)

dump(list="gompertz.makeham.par.men",file="gompertz_makeham_par_men.R")
dump(list="gompertz.makeham.est.men",file="gompertz_makeham_est_men.R")

###########
##2e: Siler, women:
###########

siler.par.wom <- list()
siler.est.wom <- list()

for(i in 1:length(countrySelect)){
	
	print(i)	

	current.Dx <- deaths.wom[[countrySelect[i]]]
	current.Nx <- exposures.wom[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=5)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.siler,lower=c(0.000001,0,0,50,0.04),
						upper=c(0.7,6,0.1,110,1), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("a","b1","c","M","b2")
	rownames(current.par) <- colnames(current.Dx)

	siler.par.wom[[countrySelect[i]]] <- current.par


	current.siler.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.siler.estimates[,year] <- siler(a1=current.par[year,1],
											b1=current.par[year,2],
											c=current.par[year,3],
											M=current.par[year,4],
											b2=current.par[year,5],
											x=current.age)
	}

	rownames(current.siler.estimates) <- current.age
	colnames(current.siler.estimates) <- rownames(current.par)

	siler.est.wom[[countrySelect[i]]] <- current.siler.estimates

}

setwd(the.dump.path)

dump(list="siler.par.wom",file="siler_par_wom.R")
dump(list="siler.est.wom",file="siler_est_wom.R")

###########
##2f: Siler, men:
###########

siler.par.men <- list()
siler.est.men <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.Dx <- deaths.men[[countrySelect[i]]]
	current.Nx <- exposures.men[[countrySelect[i]]]
	current.age <- 0:110
	current.par <- matrix(NA,nrow=ncol(current.Dx),ncol=5)

		for(time in 1:nrow(current.par)){
			current.par[time,] <- DEoptim(fn=ll.poisson.siler,lower=c(0.000001,0,0,50,0.04),
						upper=c(0.7,6,0.1,110,1), Dx=current.Dx[,time], 
						Nx=current.Nx[,time],x=current.age,
						control=DEoptim.control(trace=300))$optim$bestmem
		}

	colnames(current.par) <- c("a","b1","c","M","b2")
	rownames(current.par) <- colnames(current.Dx)

	siler.par.men[[countrySelect[i]]] <- current.par


	current.siler.estimates <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:nrow(current.par)){
		current.siler.estimates[,year] <- siler(a1=current.par[year,1],
											b1=current.par[year,2],
											c=current.par[year,3],
											M=current.par[year,4],
											b2=current.par[year,5],
											x=current.age)
	}

	rownames(current.siler.estimates) <- current.age
	colnames(current.siler.estimates) <- rownames(current.par)

	siler.est.men[[countrySelect[i]]] <- current.siler.estimates

}

setwd(the.dump.path)

dump(list="siler.par.men",file="siler_par_men.R")
dump(list="siler.est.men",file="siler_est_men.R")


###########
##2g: Gompertz, beta up constant, beta down constant, women:
###########

gompertz.est.betaUpConst.wom <- list()
gompertz.est.betaDownConst.wom <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.par <- gompertz.par.wom[[countrySelect[i]]] 
	current.years <- rownames(current.par)
	current.age <- 0:110
	current.gompertz.estimates.betaUpConst <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))
	current.gompertz.estimates.betaDownConst <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:length(current.years)){
		current.gompertz.estimates.betaUpConst[,year] <- gompertzDifferBeta(M=current.par[year,1],
														bDown=current.par[year,2],
														bUp=current.par[as.character(1985),2],
														x=current.age)
		
		current.gompertz.estimates.betaDownConst[,year] <- gompertzDifferBeta(M=current.par[year,1],
														bDown=current.par[as.character(1985),2],
														bUp=current.par[year,2],
														x=current.age)
	}

	rownames(current.gompertz.estimates.betaUpConst) <- current.age
	colnames(current.gompertz.estimates.betaUpConst) <- rownames(current.par)

	rownames(current.gompertz.estimates.betaDownConst) <- current.age
	colnames(current.gompertz.estimates.betaDownConst) <- rownames(current.par)

	gompertz.est.betaUpConst.wom[[countrySelect[i]]] <- current.gompertz.estimates.betaUpConst
	gompertz.est.betaDownConst.wom[[countrySelect[i]]] <- current.gompertz.estimates.betaDownConst

}

setwd(the.dump.path)

dump(list="gompertz.est.betaUpConst.wom",file="gompertz_est_betaUpConst_wom.R")
dump(list="gompertz.est.betaDownConst.wom",file="gompertz_est_betaDownConst_wom.R")

###########
##2h: Gompertz, beta up constant, beta down constant, men:
###########

gompertz.est.betaUpConst.men <- list()
gompertz.est.betaDownConst.men <- list()

for(i in 1:length(countrySelect)){
		
	print(i)	

	current.par <- gompertz.par.men[[countrySelect[i]]] 
	current.years <- rownames(current.par)
	current.age <- 0:110
	current.gompertz.estimates.betaUpConst <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))
	current.gompertz.estimates.betaDownConst <- matrix(NA,nrow=length(current.age),ncol=nrow(current.par))

	for(year in 1:length(current.years)){
		current.gompertz.estimates.betaUpConst[,year] <- gompertzDifferBeta(M=current.par[year,1],
														bDown=current.par[year,2],
														bUp=current.par[as.character(1985),2],
														x=current.age)
		
		current.gompertz.estimates.betaDownConst[,year] <- gompertzDifferBeta(M=current.par[year,1],
														bDown=current.par[as.character(1985),2],
														bUp=current.par[year,2],
														x=current.age)
	}

	rownames(current.gompertz.estimates.betaUpConst) <- current.age
	colnames(current.gompertz.estimates.betaUpConst) <- rownames(current.par)

	rownames(current.gompertz.estimates.betaDownConst) <- current.age
	colnames(current.gompertz.estimates.betaDownConst) <- rownames(current.par)

	gompertz.est.betaUpConst.men[[countrySelect[i]]] <- current.gompertz.estimates.betaUpConst
	gompertz.est.betaDownConst.men[[countrySelect[i]]] <- current.gompertz.estimates.betaDownConst

}

setwd(the.dump.path)

dump(list="gompertz.est.betaUpConst.men",file="gompertz_est_betaUpConst_men.R")
dump(list="gompertz.est.betaDownConst.men",file="gompertz_est_betaDownConst_men.R")

#######################################################################
#######################################################################

