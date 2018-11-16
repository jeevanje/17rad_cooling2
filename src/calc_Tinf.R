Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(gsl)
W = lambert_W0

source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source("./src/plot_params.R")
source("./src/cts_functions.R")
load("./data/band_params.Rdata")  # includes ncoarse

# params
kappamin  = min(kappa_simple)
RH        = 1
pref	  = ps
Ts        = 350
Tstrat    = 200

# Gamma calc, Coldblatt 2013, Fig. S1
Gamma   = g/Rd*log(250/350)/log(10e2/ps)
WVP0    = calc_WVP0(c(Ts,Tstrat),RH,Gamma)
Tinf    = L/Rv/log(kappamin*WVP0)

# Calc T1
Tstar   = Rd*Gamma*L/g/Rv
tau_inf = WVP0*kappa_simple*ps/pref
T1      = Tstar/lambert_W0(Tstar/Ts*tau_inf^(Rd*Gamma/g))
Tinf    = max(T1)

#============#
# romps 2016 #
#============#
Tvals = Ts:Tstrat
N     = length(Tvals)
T0	  = (Ts+Tstrat)/2
qvs   = qsat(Ts,ps)
y	  = L*qvs/Rd/T0*exp(L*qvs/Rd/T0)
f	  = L/Rv/T0^2 - Cp/Rd/T0

Tvals  = array(Tvals,dim=N)
Tvals_s= rfm_i2s(Tvals)
zvals  = Cp/g*(Ts-Tvals) + Rd*T0/g*W(y) - Rd*T0/g*W(y*exp(-f*(Ts-Tvals)))
pvals  = array(ps*exp(-g*zvals/Rd/T0),dim=N)
qsatvals= qsat(Tvals,pvals) 
qsat_alt = Rd*T0/L*W(y*exp(-f*(Ts-Tvals))) 

# Gamma calc
kstrat = which.min(Tvals-Tstrat)
zstrat = zvals[kstrat]
Gamma_romps  = -diff(Tvals)/diff(zvals)
Gamma_alt = diff(Tvals)/diff(pvals)*rfm_lni2s(pvals)/Rd/rfm_i2s(Tvals)*g
Gamma_ref = gamma_m(Tvals,pvals)

gammalim = c(0,10)  #K/km
pdf("plots/gamma_m_romps.pdf",width=5,height=5)
plot(1e3*Gamma_ref,Tvals,type="l",ylim=c(Ts,Tstrat),xlim=gammalim,
	xlab=expression(Gamma~~"(K/km)"),ylab=tabslab)
points(1e3*Gamma_romps,Tvals_s,type="l",col="red")
points(1e3*Gamma_alt,Tvals_s,type="l",col="blue")
dev.off()
