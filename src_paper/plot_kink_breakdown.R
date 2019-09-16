Rtoolsdir = "~/Dropbox/Rtools/"

library(fields)
library(ncdf4)
setwd("~/Dropbox/17rad_cooling2/src_paper/")
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("cts_functions.R")
source("plot_params.R")
load("../data_paper/band_params.Rdata")

case    = "h2o_noctm_Ts300_rh0.75_gamma7"
Gamma   = 7e-3
RH	    = 0.75
Ts	    = 300
beta    = 5.5
k_wind  = 1e5		  # window, m^-1
pref    = 5e4

#=======#
# Data  #
#=======#
ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k	    = ncvar_get(nc,"k")
p	    = ncvar_get(nc,"p")
z	    = ncvar_get(nc,"z")
tabs    = ncvar_get(nc,"tabs")
qv      = ncvar_get(nc,"q_h2o")
kref    = which.min(abs(p-pref))
dk      = k[2]-k[1]
ivec_rot = which(k<k_wind)	
nk_rot  = length(ivec_rot) 
opt     = D*ncvar_get(nc,"opt")[ivec_rot,]  			   # rot band only!
kapparef = ncvar_get(nc,"sigma_h2o",start=c(1,kref),count=c(nk_rot,1))*N_avo/m_h2o
np	    = length(p)
p_s     = rfm_lni2s(p)
qv_s    = rfm_lni2s(qv)
tabs_s  = rfm_i2s(tabs)
rho_s   = p_s/Rd/tabs_s
rhov_s  = qv_s*rho_s
pvec    = np:1
krot_rfm = k[ivec_rot]


#==========#
# Process  #
#==========#

# diagnose tp
k_tp = min(which(tabs==200))
p_tp = p[k_tp]

# WVP
wvp_rfm	= numeric(np)
wvp_rfm[np] = 0
for (k_z in (np-1):1){
	wvp_rfm[k_z] = wvp_rfm[k_z+1] + rhov_s[k_z]*diff(z)[k_z]
}
wvp_theory = calc_WVP0(tabs,RH,Gamma)*exp(-L/Rv/tabs)  # i lev

# dtrans_rfm
trans 	        = exp(-opt)
dtrans_rfm_dp2d = (trans[ ,2:np]-trans[ ,1:(np-1)])/(rep(1,times=nk_rot)%o%diff(p))  
dtrans_rfm_dp1d = apply(dtrans_rfm_dp2d,2,sum)*dk  # cm^-1/Pa

# dtrans_ssm2d
krot_ssm   = k_fit_rot
nk_rot_ssm = length(krot_ssm)
dk_rot_ssm = diff(krot_ssm)[1]
tau_ssm    = kappa_simple_h2o%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*wvp_theory)
trans_ssm  = exp(-tau_ssm)
dtrans_ssm_dp2d = (trans_ssm[1:nk_rot_ssm,2:np]-trans_ssm[1:nk_rot_ssm,1:(np-1)])/(rep(1,times=nk_rot_ssm)%o%diff(p)) #s lev
dtrans_ssm2d_dp1d = dk_rot_ssm*apply(dtrans_ssm_dp2d,2,sum)

# dtrans_ssm1d
k1_rot  = calc_k1("h2o",p,tabs,Ts,Gamma,RH)[[1]]
dtrans_ssm1d_dp1d = -l_rot*5.5/p_s
dtrans_ssm1d_dp1d[is.na(k1_rot)[-np]] <- NA


# deltak, needs lower dk
dk_ssm     = 100  # m^-1
kvals_ssm  = seq(from=k_rot,to=1e5,by=dk_ssm)
kappavals_ssm = kappa_rot*exp(-abs(kvals_ssm -k_rot)/l_rot)
tauvals_ssm  = kappavals_ssm%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*wvp_theory)
kvec	     = 1:nk_rot
deltak_rfm   = numeric(np)
deltak_ssm2d = numeric(np)
taumin = exp(-exp(1)/2)
taumax = exp(exp(1)/2)
for (k_p in 1:np){
	deltak_rfm[k_p] = dk*sum(((opt[kvec,k_p]>taumin)&(opt[kvec,k_p]<taumax))) # m^-1
	deltak_ssm2d[k_p] = dk_ssm*sum(((tauvals_ssm[ ,k_p]>taumin)&(tauvals_ssm[ ,k_p]<taumax)))
}
deltak_ssm2d[(k_tp+1):np] = NA
deltak_ssm1d = rep(l_rot*exp(1),times=np)
deltak_ssm1d[is.na(k1_rot)] = NA

# kapparef histogram
kappa_med    = median(kapparef)
lnkappa_hist = hist(log(kapparef),breaks=40,plot=FALSE)
lnkappa_mids = lnkappa_hist$mids
lnkappa_density = lnkappa_hist$density

#=======#
# Plot  #
#=======#

# Plot params
klim           = c(10,1000)
dtrans_dp2dlim = c(0,0.85e-2)  # 1/hPa
dtrans_dp1dlim = c(0,2)  	   # cm^-1/hPa
deltak_lim 	   = c(0,250)  	   # cm^-1
wvplim   	   = c(45,2e-4)    # kg/m^2
kappalim 	   = 1/wvplim      # m^2/kg
wvplab 		   = expression(D*frac(p,p[ref])*W*V*P~~"("*kg/m^2*")")

methods = c("rfm","ssm2d","ssm1d")
methodnames = c("RFM","SSM2D","SSM1D")
colvec  = c("black","red","red")
ltyvec  = c("solid","dashed","solid")
lwd     = 4
lwdvec  = c(lwd,lwd,lwd-2)
Nmethod = length(methods)

cex	   	   = 2.25
cex_leg    = 1.75
cex.legend = 1.75
pmax       = 10e4
plim       = 1e-2*c(pmax,0)
pvec       = (np-1):1
tp_col     = "gray"
tp_lty     = "dotted"

# PDF
file = "~/Dropbox/17rad_cooling2/plots_paper/kink_breakdown.pdf"
pdf(file,width=15,height=10,bg="white")
par(oma=c(0,1,0,1))
set.panel(2,3)

par(mar=c(6,5,5,5))
# dtrans_dp2d, RFM
ncoarse     = 10/(dk/1e2) 
krot 	    = coarse_grain(krot_rfm,ncoarse)
dtrans_dp2d = coarse_grain(dtrans_rfm_dp2d,ncoarse)
main = bquote("(a)"~~~-partialdiff[p]*T[tilde(nu)]~",  RFM   ("*1/hPa*")")
my.image.plot(1e-2*krot,1e-2*p_s[pvec],-1e2*dtrans_dp2d[ ,pvec],
	xlim = klim,
	ylim = plim,
	zlim = dtrans_dp2dlim,
	xlab = klab,
	ylab = plab,
	main = main,
	cex.axis = cex,
	cex.lab  = cex,
	cex.main = cex,
	cex.legend = cex.legend
) 
abline(h=1e-2*p_tp,col=tp_col,lty=tp_lty,lwd=lwd)

par(mar=c(6,8,5,1))
# dtrans_dp1d
plot(1,type="n", 
	xlim 	 = dtrans_dp1dlim, 
	ylim	 = plim,
	xlab 	 = expression(-bar(partialdiff[p]*T[tilde(nu)])~~~"("*c*m^{-1}/h*P*a*")"),
	ylab     = plab,
	main     = "(b) Integrated transmissivity gradient",
	font.main = 1,
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
for (n_method in 1:Nmethod){
	method = methods[n_method]
	dtrans_dp1d = eval(as.name(paste("dtrans",method,"dp1d",sep="_")))
	lwd   = lwdvec[n_method]
	lty   = ltyvec[n_method]
    col   = colvec[n_method]
    points(-dtrans_dp1d,1e-2*p_s,type="l",lwd=lwd,lty=lty,col=col)
} # field
legend("bottomright",methodnames,lty=ltyvec,col=colvec,lwd=lwdvec[1]-1,cex=cex_leg)
abline(h=1e-2*p_tp,col=tp_col,lty=tp_lty,lwd=lwd)

# deltak profiles
plot(1,type="n", 
	xlim 	 = deltak_lim, 
	ylim	 = plim,
	xlab 	 = expression(Delta*tilde(nu)~~~"("*c*m^{-1}*")"),
	ylab     = plab,
	main     = "(c) Effective emitting width",
	font.main = 1,
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
for (n_method in 1:Nmethod){
	method = methods[n_method]
	deltak = eval(as.name(paste("deltak",method,sep="_")))
	lwd   = lwdvec[n_method]
	lty   = ltyvec[n_method]
    col   = colvec[n_method]
    points(1e-2*deltak,1e-2*p,type="l",lwd=lwd,lty=lty,col=col)
} # field
legend("bottomleft",methodnames,lty=ltyvec,col=colvec,lwd=lwdvec[1]-1,cex=cex_leg)
abline(h=1e-2*p_tp,col=tp_col,lty=tp_lty,lwd=lwd)

# kapparef histogram
plot(lnkappa_density,exp(lnkappa_mids),type="l", lwd = 3,
	xlab 	 = "Density",
	ylab 	 = expression(kappa[ref]~~"("*m^2/kg*")"),
	ylim     = kappalim,
	log      = "y",
	main     = expression((d)~~kappa[ref]~~"distribution"),
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex 
	)

# deltak vs. WVP
plot(1,type="n", 
	xlim 	 = deltak_lim, 
	ylim	 = wvplim,
	xlab 	 = expression(Delta*tilde(nu)~~~"("*c*m^{-1}*")"),
	ylab     = wvplab,
	log      = "y",
	main     = "(e) Effective emitting width",
	font.main = 1,
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
for (n_method in 1:1){
	method = methods[n_method]
	deltak = eval(as.name(paste("deltak",method,sep="_")))
	lwd   = lwdvec[n_method]
	lty   = ltyvec[n_method]
    col   = colvec[n_method]
    points(1e-2*deltak[-np],(D*p*wvp_rfm/pref)[-np],type="l",lwd=lwd,lty=lty,col=col)
} # field
#legend("bottomleft",methodnames,lty=ltyvec,col=colvec,lwd=lwdvec[1]-1,cex=cex_leg)

dev.off()
