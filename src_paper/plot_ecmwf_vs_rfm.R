Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"gray_model.R",sep=""))
source(paste(Rtoolsdir,"/my_image_plot.R",sep=""))
source("src_paper/plot_params.R")
source("src_paper/cts_functions.R")

#=======#
# ECMWF #
#=======#
ecmwf_1dfile  = "data/ecmwf_global_mean.txt"
ecmwf_1ddata  = read.table(ecmwf_1dfile)
ecmwf_1dp     = ecmwf_1ddata[ ,1]
coo_ecmwf1d   = -ecmwf_1ddata[ ,3]

ecmwf_2dfile  = "data/ecmwf_hr_zm.nc"
ecmwf_2dnc    = nc_open(ecmwf_2dfile)
ecmwf_2dH     = ncvar_get(ecmwf_2dnc,"CSLWR")
ecmwf_2dp     = ncvar_get(ecmwf_2dnc,"dimz")
ecmwf_lat     = ncvar_get(ecmwf_2dnc,"dimy")
nlat	      = length(ecmwf_lat)

#=======#
# RFM   #
#=======#

case	= c("h2o_noctm_Ts300_rh0.75_gamma7")
ncpath  = paste("data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk      = k[2]-k[1]
p       = ncvar_get(nc,"p")
z       = ncvar_get(nc,"z")
qv      = ncvar_get(nc,"q_h2o")
tabs    = ncvar_get(nc,"tabs")
coo2d   = ncvar_get(nc,"coo")/abs(coo2d_fac)         # convert to SI 
opt     = ncvar_get(nc,"opt")
Fk      = ncvar_get(nc,"flx")					   # W/m^2/cm-1

np	    = length(p) 
p_s     = rfm_lni2s(p)
qv_s    = rfm_lni2s(qv)
tabs_s  = rfm_i2s(tabs)
coo2d_s = rfm_i2s(coo2d)
opt_s   = rfm_lni2s(opt)
cts2d_s = calc_cts2d(k,p,tabs_s,opt)
coo1d   = apply(coo2d_s,2,sum)*dk  # W/m^2/Pa
cts1d   = apply(cts2d_s,2,sum)*dk
F	    = apply(Fk,2,sum)*dk*1e-2    # W/m^2
olr	    = F[np]
Ts	    = tabs[1]
Q	    = olr-F[1]

#=======#
# Gray  #
#=======#
rho_s    = p_s/Rd/tabs_s
rhov_s   = qv_s*rho_s
wvp      = compute_wvp_k(z[-np],rhov_s,"i")  # dam i levs
kappa    = tune_kappa_Q(wvp,tabs_s,Ts,Q)
tau_gray = kappa*wvp
Bvals    = B(tabs_s)
Bs	     = B(Ts)
F_gray   = compute_Ugray(tau_gray,Bvals,Bs) - compute_Dgray(tau_gray,Bvals)
coo_gray1d = -diff(F_gray)/diff(p)
k_tau1   = which.min(abs(tau_gray-1))

#=======#
# Plot  #
#=======#

# params
p_lim       = c(1000,10) # hPa
p_fac	    = 1e-2        # SI to hPa
coo1d_lim   = c(-7,0)   # K/day
coo_ecmwf_lim = 2.75*c(-1,1) # K/day
fields1d    = c("coo","cts","coo_gray")
cex	        = 1.25
cex_leg     = 1
lwd         = 2.5
ltyvec      = c("solid","dashed","dotted")
colvec      = c("black","black","gray")
col	        = two.colors(n=64,start="blue",end="red")
lat_ticks   = c(-90,-45,-30,0,30,45,90)
p_ticks     = seq(0,1000,by=200)
# PDF
file = "plots_paper/ecmwf_vs_rfm.pdf"
pdf(file,width=10,height=5,bg="white")
par(mfrow=c(1,2),mar=c(5,6,5,3))

# ECMWF ZM
my.image.plot(sin(pi/180*ecmwf_lat[nlat:1]),ecmwf_2dp,ecmwf_2dH[nlat:1,],
#contour(ecmwf_lat[nlat:1],ecmwf_2dp,ecmwf_2dH[nlat:1,],
#	levels = seq(from=-3.5,to=1.5,by=1),
	xlab  = latlab,
	ylab  = plab,
	main  = "(a) ECMWF LW clear-sky heating rate (K/day)",
	zlim  = coo_ecmwf_lim,
	ylim  = p_lim,
	col   = col,
	cex.axis = 1e-4,
	cex.main = cex,
	cex.lab  = cex,
	cex.legend = cex_leg)
title("")
axis(side=1,at=sin(pi*lat_ticks/180),labels=lat_ticks,cex.axis=cex)
axis(side=2,at=p_ticks,labels=p_ticks,cex.axis=cex)
#mtext("K/day",side=3,line=1,outer=FALSE,adj=1,cex=1)
contour(sin(pi/180*ecmwf_lat[nlat:1]),ecmwf_2dp,ecmwf_2dH[nlat:1,],
		levels=c(-1.5,-2.5),lty="dashed",add=TRUE,lwd=2,labcex=1,col="white")

# 1D Profiles
plot(1,type="n", 
     xlim     = coo1d_lim, 
     ylim     = p_lim,
     xlab     = "H (K/day)",    
     ylab     = plab,
     main     = "(b) Heating rate profiles",
     cex.lab  = cex,
     cex.axis = cex,
     cex.main = cex)
for (n in 1:length(fields1d)){
     field = fields1d[n]
     lty   = ltyvec[n]
     var   = eval(as.name(paste(field,"1d",sep="")))
     }
points(-coo_ecmwf1d,ecmwf_1dp,type="l",lwd=lwd,lty=ltyvec[1], col=colvec[1])     
points(coo1d_fac*coo1d,p_fac*p_s,type="l",lwd=lwd,lty=ltyvec[2],col=colvec[2])  
points(coo1d_fac*coo_gray1d,p_fac*p_s,type="l",lwd=lwd,lty=ltyvec[3],col=colvec[3]) 
points(coo1d_fac*coo_gray1d[k_tau1],p_fac*p_s[k_tau1],pch=8,
	col="gray",cex=1.5) 
legend("bottomleft",legend=c("ECMWF mean","RFM 1D","Gray 1D"),lwd=lwd,
			lty=ltyvec,cex=0.9,col=colvec)

dev.off()	


