Rtoolsdir = "~/Dropbox/Rtools/"

library(ncdf4)
setwd("~/Dropbox/17rad_cooling2/src/")
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source("plot_params.R")
load("../data_paper/band_params.Rdata")

case     = "co2_noctm_Ts300_rh0.75_gamma7"
Gamma   = 7e-3
RH	    = 0.75
Ts	    = 300
D	    = 1.5

#=======#
# Data  #
#=======#

ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k	    = ncvar_get(nc,"k")
nk      = length(k)
dk	    = k[2]-k[1]
p	    = ncvar_get(nc,"p")
q	    = ncvar_get(nc,"q_co2")
np	    = length(p)
p_s     = rfm_lni2s(p)
tau_rfm = D*ncvar_get(nc,"opt")
kappa_rfm = ncvar_get(nc,"sigma_co2")*N_avo/m_co2

#==========#
# Process  #
#==========#

# kappaone
kappaone_rfm   = numeric(np)
kappaone_ssm1d = numeric(np)
taumin 		   = exp(-exp(1)/2)
taumax 		   = exp(exp(1)/2)
for (k_p in 1:np){	
	indices 		   = which((tau_rfm[ ,k_p]>taumin)&(tau_rfm[ ,k_p]<taumax))
	kappaone_rfm[k_p]  = sum(kappa_rfm[indices,k_p])/length(indices)
}
beta           = 2
kappaone_ssm1d = beta*g/q/p

#=======#
# Plot  #
#=======#

# Plot params
kappaone_lim = c(1e-1,1e3) # m^2/kg
methods     = c("rfm","ssm1d")
methodnames = c("RFM","2g/qp")
colvec  	= c("black","red")
ltyvec  	= c("solid","solid")
Nmethod 	= length(methods)

cex	   	   = 1.25
cex_leg    = 1
lwd	       = 2.5
plim       = c(1e3,1e0)
pvec       = 1:(np-2)

# PDF
file = "~/Dropbox/17rad_cooling2/plots/kappaone.pdf"
pdf(file,width=5,height=5,bg="white")
par(mfrow=c(1,1),mar=c(6,5,5,5))

# kappaone profiles
plot(1,type="n", 
	xlim 	 = kappaone_lim, 
	ylim	 = plim,
	log      = "xy",
	xlab 	 = expression(kappa[1]~~~"("*m^{2}/kg*")"),
	ylab     = plab,
	main     = "CO2 effective absorption coefficient",
	cex.lab  = cex,
	cex.axis = cex,
	cex.main = cex)
for (n_method in 1:2){
	method = methods[n_method]
	kappaone = eval(as.name(paste("kappaone",method,sep="_")))
	lty   = ltyvec[n_method]
    col   = colvec[n_method]
    points(kappaone[pvec],1e-2*p[pvec],type="l",lwd=lwd,lty=lty,col=col)
} # field
legend("topleft",methodnames,lty=ltyvec,col=colvec,lwd=lwd,cex=cex_leg)

dev.off()
