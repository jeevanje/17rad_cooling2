Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("src/plot_params.R")
load("data_paper/band_data.Rdata")
load("data_paper/band_params.Rdata")

# Data params
case_h2o  = "h2o_noctm_Ts300_rh0.75_gamma7"
case_co2  = "co2_noctm_Ts300_rh0.75_gamma7"
D         = 1.5  # diffusion parameter
RH		  = 0.75
Gamma	  = 7e-3
gases     = c("h2o","co2")
gas	      = gases[1]
case      = eval(as.name(paste("case",gas,sep="_")))

#=======#
# RFM   #
#=======#

gasname = gasnames[1]
m_mol   = eval(as.name(paste("m_",gas,sep="")))
sigmavar = paste("sigma",gas,sep="_")	
ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk	    = k[2]-k[1]
p       = ncvar_get(nc,"p")
p_s     = rfm_lni2s(p)
q       = ncvar_get(nc,"q_co2")
q_s     = rfm_i2s(q)
tabs    = ncvar_get(nc,"tabs")
tabs_s  = rfm_i2s(tabs)
Ts      = tabs[1]
Tav     = mean(range(tabs))
np      = length(p)
p_vec	= (np-1):1
nk      = length(k)
ncoarse	= 1000/dk

opt	    = D*ncvar_get(nc,"opt")
coo	    = ncvar_get(nc,"coo")/abs(coo2d_fac)         # convert to SI 
coo_s   = rfm_i2s(coo)
coo_rfm = coarse_grain(coo_s,ncoarse)
k_rfm   = coarse_grain(k,ncoarse)

#=========#
# Theory  #
#=========#

band_data_gas = band_data[[gas]]

for (j in 1:2){
   data = band_data_gas[[j]]       
   for (var in c("k_fit","kappa_fit","lk")){
       assign(paste(var,j,sep=""),data[[var]])
   }
}
k_fit     = c(k_fit1,k_fit2[-1])
dk_fit    = k_fit[2] - k_fit[1]
if (gas=="h2o"){
	WVP_s 		  = D*RH*esat(tabs_s)*Tav/L/Gamma  # Koll 2018, include D
	WVP 		  = D*RH*esat(tabs)*Tav/L/Gamma  
	kappa_theory  = c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
	                   kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr))
	tau_theory_i  = kappa_theory%o%(p/pref*WVP)                            
}
coo_theory = calc_cts2d(k_fit,p,tabs_s,tau_theory_i)  # s levs
k_theory   = k_fit

#=======#
# Plot  #
#=======#

# Plot params
# y coord
k_lim_h2o  = c(0,1450) # cm^-1  
k_lim_co2  = c(510,840) # cm^-1  
k_lim      = eval(as.name(paste("k_lim",gas,sep="_")))
p_lim      = c(1000,10) # hPa
p_fac	   = 1e-2        # SI to hPa
p_lab	   = plab
coo_lim     = c(-0.0125,0) # K/day/cm^-1
coo_units   = "(K/day/cm^-1)"
cex	       = 2
cex_leg    = 1
lwd        = 2

methods    = c("rfm","theory")
cols       = c("black","red")
ltys       = c("solid","dashed")
pvals      = c(2.5e4,5e4,7.5e4)  # Pa
Npvals     = length(pvals)
chvals     = array(dim=c(2,Npvals))

# PDF
file = "~/Dropbox/17rad_cooling2/plots_paper/chk_vs_k.pdf"
pdf(file,width=5.5,height=12,bg="white")
par(mfrow=c(3,1),mar=c(0,6,0,3),oma=c(5,0,3,0))

for (k_pval in 1:Npvals){
	pval = pvals[k_pval]
	plot(1,type="n",xlim=k_lim,ylim=coo_lim,
		 xlab = klab,
		 ylab = expression(H[tilde(nu)]*"  ("*K/day/cm^{-1}*")"),
		 main = "",
		 col.axis = "white",
		 tcl      = 0.01,
		 cex.lab  = cex)
	axis(2,cex.axis=cex,at=c(-0.01,-0.005,0))
	if (k_pval == Npvals){
		axis(1,cex.axis=cex)
		mtext(klab,1,line=3,cex=1.4)		
	} else {
		axis(1,col.axis="white")
	}
	text(1200,-0.01,paste("p=",1e-2*pval,"  hPa",sep=""),cex=cex,font=2)
	#mtext(paste("p=",1e-2*pval,"  hPa",sep=""),4,line=1,cex=1.5,font=2,las=3)
	k_p = which.min(abs(p_s-pval))
	for (k_method in 1:2){
		method    = methods[k_method]
		col       = cols[k_method]
		lty       = ltys[k_method]
		kvals     = eval(as.name(paste("k_",method,sep="")))
		coo2dvals = eval(as.name(paste("coo_",method,sep="")))
		coovals   = coo2dvals[ ,k_p]
		chval     = round(diff(k_rfm)[1]*sum(coovals)*coo1d_fac,digits=1)
		chvals[k_method,k_pval] = chval
		points(1e-2*kvals,coo2d_fac*coovals,type="l",col=col,lty=lty,lwd=lwd)
		text(1200,-2e-3-k_method*1e-3,bquote(H==.(chval)~K/d*a*y),col=col,cex=1.9)
	}
} 
dev.off()	


