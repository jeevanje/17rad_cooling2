Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("src/plot_params.R")
load("data/band_data.Rdata")
load("data/band_params.Rdata")

# Data params
case	  = "h2o_only_no_cont"
D         = 1.5  # diffusion parameter
RH		  = 0.75
Gamma	  = 7e-3
Nsec      = 3600*24   # sec/day
coo_fac   = -g/Cp*Nsec*1e2  # pppf to K/day/cm^-1 heating
lntau_fac = 1 
fac1d	  = -g/Cp*Nsec      # pppf to K/day

# Plot params
# y coord
k_lim	   = c(150,1450) # cm^-1  
p_lim      = c(1000,1) # hPa
p_lab      = plab
p_fac	   = 1e-2        # SI to hPa
z_lim      = c(1,30)     # km
z_lab      = zlab
z_fac	   = 1e-3	     # SI to km
ycoord     = "p"
for (var in c("lim","lab","fac")){	assign(paste("y_",var,sep=''),eval(as.name(paste(ycoord,"_",var,sep=''))))
}

# fields
fields 	    = c("kappa","coo","lntau")
coo_lim     = c(-0.014,0.014) # K/day/cm^-1
H_lim       = c(-3,0)   #K/day
lntau_lim   = c(-20,20)
coo_units   = "(K/day/cm^-1)"
lntau_units = ""
kappa_lab   = expression(kappa[ref]~~"("*m^2/kg*")")
kappa_lim   = c(1e-4,1e4) 
cex	       = 2.5
cex_leg    = 2.25
lwd        = 2
ltyvec     = c("solid","dashed")
colvec     = two.colors(n=64,start="blue",end="red",middle="white")
nrun   	   = 10  # every 10 cm^-1
ncoarse	   = 10

# Functions
plot_kappa = function(k,kappa,main){
	plot(1e-2*k,kappa,type="l",
		xlim 	 = k_lim, 
		ylim	 = kappa_lim,
		xlab 	 = klab,
		ylab     = kappa_lab,
		log      = "y",
		main	 = main,
		lwd		 = lwd,
		cex.lab  = cex,
		cex.axis = cex,
		cex.main = cex)
}


# PDF
gas     = substr(case,1,3)	
file = paste("~/Dropbox/17rad_cooling2/plots/",gas,"_rfm_theory.pdf",
			sep="")
pdf(file,width=17,height=9,bg="white")
par(mfrow=c(2,3),mar=c(5,6,5,12))

#=======#
# RFM   #
#=======#

gasname = gasnames[1]
m_mol   = eval(as.name(paste("m_",gas,sep="")))	
ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk	    = k[2]-k[1]
p       = ncvar_get(nc,"p")
p_s     = rfm_lni2s(p)
z       = ncvar_get(nc,"z")
y  	    = eval(as.name(ycoord))
y_s     = rfm_i2s(y)
tabs    = ncvar_get(nc,"tabs")
tabs_s  = rfm_i2s(tabs)
Ts      = tabs[1]
np      = length(p)
p_vec	= (np-1):1
z_vec   = 1:(np-1)
y_vec   = eval(as.name(paste(ycoord,"_vec",sep="")))
nk      = length(k)
k_coarse= coarse_grain(k,ncoarse)
nk_coarse =length(k_coarse)
dk_coarse = k_coarse[2]-k_coarse[1]

kappa   = ncvar_get(nc,"sigma",start=c(1,1),count=c(nk,1))*N_avo/m_mol # m^2/kg
opt     = D*ncvar_get(nc,"opt")  # plot D*tau since that's what's relevant
opt_s   = rfm_lni2s(opt)
lntau   = log(opt)
lntau_s = log(opt_s)
coo	    = ncvar_get(nc,"coo")/abs(coo_fac)         # convert to SI 
coo_s   = rfm_i2s(coo)
lntau_rfm= coarse_grain(lntau_s,ncoarse)
p1_rfm  = array(dim=c(nk))
for (m in 1:length(p1_rfm)){
    p1_rfm[m] = p_s[which.min(abs(lntau_s[m,]))]
}
p1_coarse = coarse_grain(p1_rfm,ncoarse)

# Kappa
kappa_coarse = coarse_grain(kappa,ncoarse)
main 		 = "(a)  absorption spectrum, RFM"
plot_kappa(k_coarse,kappa_coarse,main)

# lntau
main		 = bquote("(b)"~~ln~D*tau[k]*", RFM")
plot_2dfield(k_coarse,lntau_rfm,1,lntau_lim,main)
points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=1)

# H_k
coo_coarse    = coarse_grain(coo_s,ncoarse)
main		  = bquote("(c)"~~H[k]*", RFM  ("~K/day/cm^{-1}~")")
plot_2dfield(k_coarse,coo_coarse,coo2d_fac,coo2d_lim,main)
#points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=1)

#H
#H_rfm = fac1d*apply(coo_coarse,2,sum)*dk_coarse
#plot(H_rfm,1e-2*p_s,type="l",xlim=H_lim,ylim=p_lim)

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
k_fit     = c(k_fit1,k_fit2)
dk_fit    = k_fit[2] - k_fit[1]
kappa_fit = c(kappa_fit1,kappa_fit2)
main 	  = "(d)  absorption spectrum, fit"
plot_kappa(k_fit,kappa_fit,main)   

# lntau
WVP_s 		  = D*RH*esat(tabs_s)*tabs_s/L/Gamma  # Koll 2018, include D
WVP 		  = D*RH*esat(tabs)*tabs/L/Gamma  
k_fit         = c(k_fit_rot,k_fit_vr[-1])
nk_fit        = length(k_fit)
kappa_theory  = c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
                   kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr))
tau_theory_i  = kappa_theory%o%(p/pref*WVP)                            
tau_theory_s  = kappa_theory%o%(p_s/pref*WVP_s)         
lntau_theory  = log(tau_theory_s)
k1_rot  = k_rot + l_rot*(log(kappa_rot*WVP_s)+ log(p_s/pref))
k1_vr   = k_vr  - l_vr*(log(kappa_vr*WVP_s)  + log(p_s/pref))
k1_rot[k1_rot<k_rot] <- NA
k1_vr[k1_vr>k_vr]    <- NA
main = bquote("(e)"~~ln~D*tau[k]*", Theory")
plot_2dfield(k_fit,lntau_theory,1,lntau_lim,main)
for (band in c("rot","vr")){
     k1  = eval(as.name(paste("k1_",band,sep="")))
     points(1e-2*k1,1e-2*p_s,type="l",lty="dashed",lwd=lwd)
}

# CTS  (s lev)
dtaudp  = (tau_theory_i[ ,-1]-tau_theory_i[ ,-np])/(rep(1,times=nk_fit)%o%diff(p)) #s lev
B       = outer(k_fit,tabs_s,planck_k)
trans   = exp(-tau_theory_s)
cts2d_theory = pi*B*trans*dtaudp
main		 = bquote("(f)"~~H[k]*", Theory  ("~K/day/cm^{-1}~")")
plot_2dfield(k_fit,cts2d_theory,coo2d_fac,coo2d_lim,main)
for (band in c("rot","vr")){
     k1  = eval(as.name(paste("k1_",band,sep="")))
     #points(1e-2*k1,1e-2*p_s,type="l",lty="dashed",lwd=lwd)
}

#H_theory
#H_theory = fac1d*apply(cts2d_theory,2,sum)*dk_fit
#plot(H_theory,1e-2*p_s,type="l",xlim=H_lim,ylim=p_lim)
 
dev.off()	


