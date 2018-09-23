Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source("./src/plot_params.R")
load("./data/band_data.Rdata")  # includes ncoarse
load("./data/band_params.Rdata")  # includes ncoarse

# Data params
cases	  = c("h2o_only_no_cont","co2_only_simple_atm")
Ncase     = length(cases)
betas	  = c(5.7,2.7)
D         = 1.5  # diffusion parameter
Nsec      = 3600*24   # sec/day
coo_fac   = -g/Cp*Nsec*1e2  # pppf to K/day/cm^-1 heating
cts_fac   = coo_fac         # pppf to K/day/cm^-1 heating
lntau_fac = 1 
fac1d	  = -g/Cp*Nsec      # pppf to K/day
RH	      = 0.75
Gamma	  = 7e-3            # K/km

# Plot params
# y coord
p_lim      = c(1000,1) # hPa
p_lab      = "p (hPa)"
p_fac	   = 1e-2        # SI to hPa
z_lim      = c(1,30)     # km
z_lab      = "z (km)"
z_fac	   = 1e-3	     # SI to km
ycoord     = "p"
for (var in c("lim","lab","fac")){
	assign(paste("y_",var,sep=''),eval(as.name(paste(ycoord,"_",var,sep=''))))
}

# fields
fields1d    = c("cts","cts_simple")
coo1d_lims  = list(c(-4,0),c(-1,0.1))   # K/day
cex	       = 1.5
cex_leg    = 2
lwd        = 2
ltyvec     = c("solid","dashed")

# PDF
file = "~/Dropbox/17rad_cooling2/plots/cts_theory.pdf"
pdf(file,width=10,height=5,bg="white")
par(mfrow=c(1,2),mar=c(5,5,5,3))

for (i in 1:Ncase){
	#=======#
	# Data  #
	#=======#

	case    = cases[i]
	gas     = substr(case,1,3)
	gasname = gasnames[i]	
	beta    = betas[i]
	m_mol   = eval(as.name(paste("m_",gas,sep="")))	
	ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk      = k[2]-k[1]
	p       = ncvar_get(nc,"p")
	np      = length(p)
	z       = ncvar_get(nc,"z")
	y  	    = eval(as.name(ycoord))
	tabs    = ncvar_get(nc,"tabs")
	Ts	    = tabs[1]
	Tstrat  = tabs[np]
	Tav     = 0.5*(Ts+Tstrat)
	q       = ncvar_get(nc,paste("q_",gas,sep=""))
	rho_gas = p/Rd/tabs*q
	np_s    = np - 1
	p_vec	= np_s:1
	z_vec   = 1:np_s
	y_vec   = eval(as.name(paste(ycoord,"_vec",sep="")))
	nk      = length(k)
	m_strat = min(which(tabs == 200))
	p_strat = p[m_strat]
	coo2d   = ncvar_get(nc,"coo")/abs(coo_fac)  # convert to SI 
	opt     = D*ncvar_get(nc,"opt")
	flx2d   = ncvar_get(nc,"flx")*1e-2          # W/m^2/m^-1
	
	#=========#
	# Process #
	#=========#
	
	# get cts1d
	p_s      = rfm_lni2s(p)
	q_s      = rfm_lni2s(q)
	y_s      = rfm_i2s(y)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	opt_s    = rfm_lni2s(opt)
	lntau2d_s= log(opt_s)
	dtaudp   = (opt[ ,2:np]-opt[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	B_s      = outer(k,tabs_s,planck_k)
	
	cts2d_s  = pi*B_s*exp(-opt_s)*dtaudp  #s, p_s, W/m^2/Pa/m^-1, pppf units
	cts1d    = apply(cts2d_s,2,sum)*dk
	
	# construct cts_simple
	cts_simple_bands = array(dim=c(np_s,2))
	band_data_gas    = band_data[[gas]]
	if (gas == "h2o"){
		WVP0 = RH*einf*Tav/L/Gamma
		k1_rot  = k_rot + l_rot*(log(kappa_rot*WVP0)+log(ps/pref)+
								 g/Rd/Gamma*log(tabs_s/Ts) - L/Rv/tabs_s)  # s levs
		k1_vr   = k_vr  - l_vr*(log(kappa_vr*WVP0)+log(ps/pref)+
								 g/Rd/Gamma*log(tabs_s/Ts) - L/Rv/tabs_s)
		k1_rot[k1_rot<k_rot] <- NA
		k1_vr[k1_vr>k_vr]    <- NA
		cts_simple_bands[ ,1] = pi*planck_k(k1_rot,tabs_s)*l_rot*beta/p_s
		cts_simple_bands[ ,2] = pi*planck_k(k1_vr,tabs_s)*l_vr*beta/p_s
	} else if (gas=="co2"){
		k1_P = k_Q - l_Q*log(q_s*p_s^2*kappa_Q/2/g/pref)
		k1_R = k_Q + l_Q*log(q_s*p_s^2*kappa_Q/2/g/pref)
		k1_P[k1_P>k_Q] <- NA
		k1_R[k1_R<k_Q] <- NA
		cts_simple_bands[ ,1] = pi*planck_k(k1_P,tabs_s)*l_Q*beta/p_s
		cts_simple_bands[ ,2] = pi*planck_k(k1_R,tabs_s)*l_Q*beta/p_s
	}
	cts_simple1d = apply(cts_simple_bands,1,sum,na.rm=TRUE)
	cts_simple1d[cts_simple1d==0] <- NA
	
	#=======#
	# Plot  #
	#=======#

	coo1d_lim = coo1d_lims[[i]]
	m_tp      = min(which(tabs_s == 200))

	# 1Dfields 
    plot(1,type="n", 
        xlim     = coo1d_lim, 
        ylim     = y_lim,
        xlab     = "H (K/day)",    
        ylab     = y_lab,
        main     = paste(gasname," cooling-to-space",sep=""),
        cex.lab  = cex,
        cex.axis = cex,
        cex.main = cex)
    for (n in 1:length(fields1d)){
        field = fields1d[n]
        lty   = ltyvec[n]
        var   = eval(as.name(paste(field,"1d",sep="")))
        points(fac1d*var,y_fac*y_s,type="l",lwd=2,lty=lty)  
     }
    #points(cts_simple_bands[ ,1]*fac1d,y_fac*y_s,
    #		type="l",lwd=2,lty=lty,col="gray")  
	abline(v=0,lty="dashed",lwd=2,col="gray")
	abline(h=1e-2*p_strat,lty="dotted",lwd=2,col="gray")
	legend("bottomleft",c("RFM","Theory"),lty=ltyvec,lwd=2,cex=1.25)
}    # cases

dev.off()	

# estimate H2O cooling:
beta  = betas[1]
kest  = 12  # p=500 hPa, T=260
Test  = tabs[kest]
pest  = p[kest]
B_rot = planck_k(k1_rot[kest],Test)
B_vr  = planck_k(k1_vr[kest],Test)
H_h2o = fac1d*(pi*B_rot*l_rot + pi*B_vr*l_vr)*beta/pest