Rtoolsdir = "~/Dropbox/Rtools/"

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("~/Dropbox/17rad_cooling2/src/plot_params.R")  
load("~/Dropbox/17rad_cooling2/data/band_data.Rdata")  # includes ncoarse
load("~/Dropbox/17rad_cooling2/data/band_params.Rdata")  # includes ncoarse

# Data params
cases     = c("h2o_only_no_cont","co2_only_simple_atm")
Ncase     = length(cases)
gas_bands = list(c("rot","vr"),c("P","R"))
Gamma	  = 7e-3  # K/m
RH	      = 0.75

# Plot params
klim_h2o  = c(150,1450)  # cm^-1
klim_co2  = c(510,840)
plim      = c(1000,1) # hPa
pfac      = 1e-2        # SI to hPa
lntaulim  = c(-20,20)
cex	      = 2
cex_leg   = 1.5
lwd       = 1.5
#ncoarse   = 10
panels    = c("(a)","(b)","(c)","(d)")

plot_tau = function(kvals,pvals,lntau,main,klim){
        my.image.plot(1e-2*kvals,1e-2*pvals[p_vec],lntau[ ,p_vec],
                        xlim = klim,
                        ylim = plim,
                        zlim = lntaulim,
                        xlab = klab,
                        ylab = plab,
                        main = main,
                        cex.lab  = cex,
	                    cex.axis = cex,
                        cex.main = cex,
                        cex.legend = cex_leg
        )
}

# PDF
file = "~/Dropbox/17rad_cooling2/plots/tau_rfm_theory.pdf"
pdf(file,width=12,height=10,bg="white")
par(mfrow=c(2,2),mar=c(5,7,5,9))

for (i in 1:2){
        case  = cases[i]

        #=======#
        # Data  #
        #=======#

        gas     = substr(case,1,3)      
 		gasname = gasnames[i]
        m_mol   = eval(as.name(paste("m_",gas,sep=""))) 
        ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
        nc      = nc_open(ncpath)
        k_rfm   = ncvar_get(nc,"k")
        p       = ncvar_get(nc,"p")
        np      = length(p)
        tabs    = ncvar_get(nc,"tabs")
		Ts	    = tabs[1]
		Tstrat  = tabs[np]
		Tav     = 0.5*(Ts+Tstrat)
        q       = ncvar_get(nc,paste("q_",gas,sep=""))
        p_vec   = (np-1):1
		lntau   = log(ncvar_get(nc,"opt")[,1:(np-1)])
        lntau_rfm = coarse_grain(lntau,ncoarse)   # no D. Coarse grain after ln!

        k_coarse   = coarse_grain(k_rfm,ncoarse)
		nk_coarse  = length(k_coarse)
        p1_rfm     = array(dim=c(nk_coarse))  
	    for (m in 1:nk_coarse){
	        p1_rfm[m] = p[which.min(abs(lntau_rfm[m,]))]
	    }
        
		#========#
		# Theory #
		#========#
		
		if (gas == "h2o"){			
			WVP0	  = RH*einf*Tav/L/Gamma
			k_fit   	  = c(k_fit_rot,k_fit_vr[-1])
			nk_fit	      = length(k_fit)
			kappa_theory  = c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
							  kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr))
 			
 			lntau_theory  = log(kappa_theory%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*WVP0*exp(-L/Rv/tabs))) 			
			k1_rot  = k_rot + l_rot*(log(kappa_rot*WVP0)+log(ps/pref)+
									 g/Rd/Gamma*log(tabs/Ts) - L/Rv/tabs)
			k1_vr   = k_vr  - l_vr*(log(kappa_vr*WVP0)+log(ps/pref)+
									 g/Rd/Gamma*log(tabs/Ts) - L/Rv/tabs)
			k1_rot[k1_rot<k_rot] <- NA
			k1_vr[k1_vr>k_vr]    <- NA
		} else if (gas == "co2"){
			k_fit   	  = c(k_fit_P,k_fit_R[-1])
			nk_fit	      = length(k_fit)
			kappa_theory  = kappa_Q*exp(-abs(k_fit-k_Q)/l_Q) 			
 			lntau_theory  = log(0.5*kappa_theory%o%(q*p^2/(g*pref)))
			k1_P = k_Q - l_Q*log(q*p^2*kappa_Q/2/g/pref)
			k1_R = k_Q + l_Q*log(q*p^2*kappa_Q/2/g/pref)
			k1_P[k1_P>k_Q] <- NA
			k1_R[k1_R<k_Q] <- NA
		}

        p1_theory     = array(dim=c(nk_fit))  
	    for (m in 1:nk_fit){
	        p1_theory[m] = p[which.min(abs(lntau_theory[m,]))]
	    }

        #=======#
        # Plot  #
        #=======#
	
		klim  = eval(as.name(paste("klim_",gas,sep="")))
		bands = gas_bands[[i]]

		# RFM
        plot_tau(k_coarse,p,lntau_rfm,
        	     bquote(.(panels[2*i-1])~~.(gasname)~~ln~tau[k]~", RFM"),klim)     
	    points(1e-2*k_coarse,1e-2*p1_rfm,type="l",lty="dashed",lwd=lwd)

		# Theory
        plot_tau(k_fit,p,lntau_theory,
        		bquote(.(panels[2*i])~~.(gasname)~~ln~tau[k]~", theory"),
        		klim)     
		for (band in bands){
			k1  = eval(as.name(paste("k1_",band,sep="")))
			points(1e-2*k1,1e-2*p,type="l",lty="dashed",lwd=lwd)
		}
}
dev.off()
