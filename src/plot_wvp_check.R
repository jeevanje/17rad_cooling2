library(ncdf)

source("~/Dropbox/Rtools/plot_tools.R")
source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/calculus_tools.R")

Gamma  = 7e-3  # K/km
RH	   = 0.75
ncpath = "~/Dropbox/17rad_cooling2/data/h2o_only_no_cont.nc"
nc	   = open.ncdf(ncpath)
z      = get.var.ncdf(nc,"z")
zint   = zinterp(z)
nz	   = length(z)
tabs   = get.var.ncdf(nc,"tabs")
Ts     = tabs[1]
Tstrat = min(tabs)
Tav    = 0.5*(Ts+Tstrat)
p	   = get.var.ncdf(nc,"p")
rho    = p/Rd/tabs
qv	   = get.var.ncdf(nc,"q_h2o")
dzvec  = c(diff(zint),zint[nz]-zint[nz-1])
rhov   = qv*rho
#kmin   = which.min(tabs)
zvec   = 1:nz

# Compute WVP
wvp_rfm	 = numeric(nz)	 				#ssi
wvp_est	 = numeric(nz)	 				#ssi
wvp_rfm[nz] = rhov[nz]*dzvec[nz]
for (k in (nz-1):1){
		wvp_rfm[k] = wvp_rfm[k+1] + dzvec[k]*rhov[k]
} 
wvp0    = RH*einf*Tav/L/Gamma
wvp_est = wvp0*exp(-L/Rv/tabs)

logvec = c("","x")
posvec = c("topright","bottomleft")
cex    = 1.25

# Begin plot
pdf("~/Dropbox/17rad_cooling2/plots/wvp_check.pdf",width=9,height=5)
par(mfrow=c(1,2),mar=c(5,5,5,3))
for (i in 1:2){
	log = logvec[i]
	pos = posvec[i]
	plot(wvp_rfm[zvec],tabs[zvec],
		  ylim = rev(range(tabs)),
		  xlab = expression("WVP  ("~kg~m^-2~")"),
		  ylab = "Temperature (K)",
		  main = "Water vapor path",
		  type = "l",
		  log  = log,
		  lwd  = 2,
		  cex.lab  = cex,
		  cex.axis = cex,
		  cex.main = cex,
		  )
	points(wvp_est[zvec],tabs[zvec],type="l",lwd=2,lty="dashed")
	legend(pos,legend=c("RFM","Theory"),
			lty=c("solid","dashed"),lwd=2,cex=1.25)
}
dev.off()
