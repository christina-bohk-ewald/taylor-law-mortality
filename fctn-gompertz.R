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

gompertz <- function(M,b,x){
	p1 <- b*exp(b*(x-M))
	return(p1)
}

ll.poisson.gompertz <- function(theta,Dx,Nx,x){
	M <- theta[1]
	b <- theta[2]
	out <- -sum(Dx*log(gompertz(x=x,M=M,b=b))-gompertz(x=x,M=M,b=b)*Nx)
	return(out)
}

gompertzDifferBeta <- function(M,bDown,bUp,x){
	p1 <- bDown*exp(bUp*(x-M))
	return(p1)
}