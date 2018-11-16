Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(ncdf4)
library(gsl)

source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source("./src/plot_params.R")
source("./src/cts_functions.R")
load("./data/band_params.Rdata")  # includes ncoarse

# params
case      = "h2o_only_no_cont"
RH        = 0.75
Gamma     = 7e-3            # K/km

# base data
ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k") 
p       = ncvar_get(nc,"p") 
Fk      = ncvar_get(nc,"flx")  # net upward flux, W/m^2/cm^-1
tabs    = ncvar_get(nc,"tabs")

# cont data
ncpath_cont  = "~/Dropbox/17rad_cooling2/data/h2o_only_simple_atm.nc"
nc_cont      = nc_open(ncpath_cont)
Fk_cont      = ncvar_get(nc_cont,"flx")  # net upward flux, W/m^2/cm^-1

# process
np		= length(tabs)
nk		= length(k)
WVP0    = calc_WVP0(tabs,RH,Gamma)
OLRk    = array(Fk[ ,np],dim=c(nk))
OLRk_cont= array(Fk_cont[ ,np],dim=c(nk))
Ts      = tabs[1]
Tstar   = Rd*Gamma*L/g/Rv
tau_inf = WVP0*kappa_simple*ps/pref
T1      = Tstar/lambert_W0(Tstar/Ts*tau_inf^(Rd*Gamma/g))
kwindvals= k_fit[c(min(which(T1>Ts)),max(which(T1>Ts)))]
T1[T1>Ts] <- Ts

#validate
# k1_list = calc_k1("h2o",p,tabs,Ts,Gamma,RH)
# k1_rot  = k1_list[[1]]
# k1_vr   = k1_list[[2]]
# plot(1,type="n",xlim=1e-2*range(k),ylim=rev(range(tabs)),xlab=klab,ylab=tabslab)
# points(1e-2*k_fit,T1,type="l",col="red")
# points(1e-2*k1_rot,tabs,type="l",col="blue")
# points(1e-2*k1_vr,tabs,type="l",col="blue") # works !

# plot prep
OLRk_coarse = coarse_grain(OLRk,ncoarse)
OLRk_cont_coarse = coarse_grain(OLRk_cont,ncoarse)
k_coarse    = coarse_grain(k,ncoarse)
cex 		= 1.25
lwd 		= 2
ltyvec      = c("solid","dotted","solid")
colvec      = c("black","black","red")
windlaby    = 0.01
windcol     = "darkgray"
windlty     = "dashed"

pdf("./plots/olr.pdf",width=6,height=5)
par(mar=c(5,5,4,2))
plot(1,type="n",xlim=1e-2*range(k),ylim=range(OLRk),xlab=klab,ylab=OLRklab,
	main     = "Spectrally-resolved OLR",
	cex.main = cex,
	cex.axis = cex,
	cex.lab  = cex)
# base	
points(1e-2*k_coarse,OLRk_coarse,type="l",col=colvec[1],
	  lty=ltyvec[1],lwd=lwd)
# w/cont
points(1e-2*k_coarse,OLRk_cont_coarse,type="l",col=colvec[2],
	  lty=ltyvec[2],lwd=lwd-1)	  
# Theory
points(1e-2*k_fit,1e2*pi*planck_k(k_fit,T1),type="l",col=colvec[3],
	  lty=ltyvec[3],lwd=lwd+1)
#points(1e-2*k,1e2*pi*planck_k(k,300),type="l",col="black",lty="dashed",lwd=2)
#points(1e-2*k,1e2*pi*planck_k(k,200),type="l",col="blue",lty="dashed",lwd=2)
#abline(v=1e-2*kwindvals,col="gray",lwd=1.5)
# Window
lines(rep(1e-2*kwindvals[1],times=2),y=c(-1,pi*1e2*planck_k(kwindvals[1],Ts)),
	col=windcol,lwd=1.5,lty=windlty)
lines(rep(1e-2*kwindvals[2],times=2),y=c(-1,pi*1e2*planck_k(kwindvals[2],Ts)),
	col=windcol,lwd=1.5,lty=windlty)
arrows(1e-2*kwindvals[1],windlaby,1e-2*kwindvals[2],code=3,
		length=0.1,col=windcol)
text(1e-2*mean(kwindvals),windlaby+0.02,"window",col=windcol)
legend("topleft",c("Base","w/cont","Theory"),
	lty=ltyvec,col=colvec,lwd=lwd,cex=1)
dev.off()
