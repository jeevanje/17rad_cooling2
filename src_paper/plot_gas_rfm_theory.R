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

# Plot params
# y coord
fields 	   = c("kappa","coo","lntau")
k_lim_h2o  = c(150,1450) # cm^-1  
k_lim_co2  = c(510,840) # cm^-1  
k_lim      = eval(as.name(paste("k_lim",gas,sep="_")))
p_lim      = c(1000,10) # hPa
p_fac	   = 1e-2        # SI to hPa
p_lab	   = plab
log_h2o    = ""
log_co2    = "y"
log        = eval(as.name(paste("log",gas,sep="_")))
coo_lim     = c(-0.014,0.014) # K/day/cm^-1
#coo2d_lim     = c(-5e1,0) # K/day/cm^-1
lntau_lim   = c(-20,20)
coo_units   = "(K/day/cm^-1)"
lntau_units = ""
kappa_lab   = expression(kappa[ref]~~"("*m^2/kg*")")
kappa_lim   = c(1e-4,1e3) 
cex	       = 2.5
cex_leg    = 2.25
lwd        = 2
colvec     = two.colors(n=64,start="blue",end="red",middle="white")
taumax     = 1

# Functions
plot_kappa = function(k,kappa,main){
	plot(1e-2*k,kappa,type="l",
		xlim 	 = k_lim, 
		ylim	 = kappa_lim,
		xlab 	 = "",
		ylab     = kappa_lab,
		log      = "y",
		main	 = "",
		lwd		 = lwd,
		#lab      = c(5,2,5),
		#yaxp     = c(1e-3,1e3,1),
		cex.lab  = cex,
		cex.axis = cex,
		cex.main = cex)
	mtext(klab,side=1,line=4,outer=FALSE,cex=1.75)	
	mtext(main,side=3,line=2,outer=FALSE,cex=1.75)	
}


# PDF
file = paste("~/Dropbox/17rad_cooling2/plots/",gas,"_rfm_theory.pdf",
			sep="")
pdf(file,width=17,height=9,bg="white")
par(mfrow=c(2,3),mar=c(5.5,6,5,12))

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
np      = length(p)
p_vec	= (np-1):1
nk      = length(k)
ncoarse	= 1000/dk
k_coarse= coarse_grain(k,ncoarse)
nk_coarse =length(k_coarse)
dk_coarse = k_coarse[2]-k_coarse[1]

kappa   = ncvar_get(nc,sigmavar,start=c(1,1),count=c(nk,1))*N_avo/m_mol # m^2/kg
opt     = D*ncvar_get(nc,"opt")  # plot D*tau since that's what's relevant
opt_s   = rfm_lni2s(opt)
lntau   = log(opt)
lntau_s = log(opt_s)
coo	    = ncvar_get(nc,"coo")/abs(coo2d_fac)         # convert to SI 
coo_s   = rfm_i2s(coo)
lntau_rfm= coarse_grain(lntau_s,ncoarse)
p1_rfm  = array(dim=c(nk))
for (m in 1:length(p1_rfm)){
    p1_rfm[m] = p_s[which.min(abs(lntau_s[m,]-log(taumax)))]
}
p1_coarse = coarse_grain(p1_rfm,ncoarse)

# Kappa
kappa_coarse = exp(coarse_grain(log(kappa),ncoarse))
main 		 = "(a)  absorption spectrum, RFM"
plot_kappa(k_coarse,kappa_coarse,main)

# lntau
main		 = bquote("(b)"~~ln~D*tau[k]*", RFM")
plot_2dfield(k_coarse,lntau_rfm,1,lntau_lim,main,log=log)
points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=1)

# H_k
coo_coarse    = coarse_grain(coo_s,ncoarse)
main		  = bquote("(c)"~~H[k]*", RFM  ("*K/day/cm^{-1}*")")
plot_2dfield(k_coarse,coo_coarse,coo2d_fac,coo2d_lim,main,log=log)
#points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=1)


#=========#
# Theory  #
#=========#

band_data_gas = band_data[[gas]]

# kappa
for (j in 1:2){
   data = band_data_gas[[j]]       
   for (var in c("k_fit","kappa_fit","lk")){
       assign(paste(var,j,sep=""),data[[var]])
   }
}
k_fit     = c(k_fit1,k_fit2[-1])
dk_fit    = k_fit[2] - k_fit[1]
kappa_fit = c(kappa_fit1,kappa_fit2[-1])
main 	  = "(d)  absorption spectrum, SSM2D"
plot_kappa(k_fit,kappa_fit,main)   
if (gas=="h2o"){
	text(500,10,"rot",cex=cex)
	text(1200,1e0,"v-r",cex=cex)
}
# lntau
if (gas=="h2o"){
	WVP_s 		  = D*RH*esat(tabs_s)*tabs_s/L/Gamma  # Koll 2018, include D
	WVP 		  = D*RH*esat(tabs)*tabs/L/Gamma  
	kappa_theory  = c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
	                   kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr))
	tau_theory_i  = kappa_theory%o%(p/pref*WVP)                            
	tau_theory_s  = kappa_theory%o%(p_s/pref*WVP_s)         
	k1_1  		  = k_rot + l_rot*(log(kappa_rot*WVP_s)+ log(p_s/pref))
	k1_2  		  = k_vr  - l_vr*(log(kappa_vr*WVP_s)  + log(p_s/pref))
	k1_1[k1_1<k_rot] <- NA
	k1_2[k1_2>k_vr]   <- NA
} else if (gas=="co2"){
	kappa_theory  = kappa_Q*exp(-abs(k_fit-k_Q)/l_Q)         
	tau_theory_i  = 0.5*kappa_theory%o%(q*p^2/(g*pref))
	tau_theory_s  = 0.5*kappa_theory%o%(q_s*p_s^2/(g*pref))
	k1_1 		  = k_Q - l_Q*log(q_s*p_s^2*kappa_Q/2/g/pref)
	k1_2 		  = k_Q + l_Q*log(q_s*p_s^2*kappa_Q/2/g/pref)
	k1_1[k1_1>k_Q] <- NA
	k1_2[k1_2<k_Q] <- NA
}
lntau_theory  = log(tau_theory_s)
main = bquote("(e)"~~ln~D*tau[k]*", SSM2D")
plot_2dfield(k_fit,lntau_theory,1,lntau_lim,main,log=log)
for (j in 1:2){
     k1  = eval(as.name(paste("k1_",j,sep="")))
     points(1e-2*k1,1e-2*p_s,type="l",lty="dashed",lwd=lwd)
}

# CTS  (s lev)
nk_fit  = length(k_fit)
dtaudp  = (tau_theory_i[ ,-1]-tau_theory_i[ ,-np])/(rep(1,times=nk_fit)%o%diff(p)) #s lev
B       = outer(k_fit,tabs_s,planck_k)
#B       = rep(1,times=nk_fit)%o%rep(1,times=np-1)
trans   = exp(-tau_theory_s)
cts2d_theory = pi*B*trans*dtaudp
main		 = bquote("(f)"~~H[k]*", SSM2D  ("*K/day/cm^{-1}*")")
plot_2dfield(k_fit,cts2d_theory,coo2d_fac,coo2d_lim,main,log=log)
for (j in 1:2){
     k1  = eval(as.name(paste("k1_",j,sep="")))
     #points(1e-2*k1,1e-2*p_s,type="l",lty="dashed",lwd=lwd)
}
 
dev.off()	


