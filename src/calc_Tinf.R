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

# Gamma calc, Goldblatt 2013, Fig. S1
Gamma_gold   = g/Rd*log(250/350)/log(10e2/ps)
# Romps, from plot_gamma_m_romps.R
Gamma_romps	 = 2e-3
Gamma   = Gamma_romps
WVP0    = calc_WVP0(c(Ts,Tstrat),RH,Gamma)
Tinf    = L/Rv/log(kappamin*WVP0)

# Calc T1, old
# Tstar   = Rd*Gamma*L/g/Rv
# tau_inf = WVP0*kappa_simple*ps/pref
# T1      = Tstar/lambert_W0(Tstar/Ts*tau_inf^(Rd*Gamma/g))
# Tinf    = max(T1)

