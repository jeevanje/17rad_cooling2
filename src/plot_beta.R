library(fields)
library(ncdf4)
source("~/Dropbox/Rtools/rfm_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")

cases     = c("co2_only_isotherm_atm","co2_only_simple_atm",
			  "h2o_no_cont_isotherm_atm","h2o_only_no_cont")
panels    = c("a","b","c","d")      	  
casenames = c("CO2, no T-scaling","CO2","H2O, no T-scaling","H2O")

#======#
# Data #
#======#

for (case in cases){
    ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
    nc	    = nc_open(ncpath)
    opt     = ncvar_get(nc,"opt")
    k       = ncvar_get(nc,"k")
    nk	    = length(k)
    z       = ncvar_get(nc,"z")
    nz      = length(z)
    p       = ncvar_get(nc,"p")
    p_s     = rfm_lni2s(p)
    opt_s   = rfm_lni2s(opt)
    beta    = (log(opt_s[ ,2:(nz-1)])-log(opt_s[ ,1:(nz-2)]))/(rep(1,times=nk)%o%diff(log(p_s)))
    z1 	    = array(dim=c(nk))
    beta1   = array(dim=c(nk))
    for (i in 1:nk){
    	m_i      = which.min(abs(opt[i,]-1))
    	z1[i]    = z[m_i]
    	m_i	     = max(1,m_i-1)
    	m_i  	 = min(m_i,nz-2)
		beta1[i] = beta[i,m_i]
    }
    assign(paste("beta_",case,sep=""),beta)
    assign(paste("beta1_",case,sep=""),beta1)
    assign(paste("k_",case,sep=""),k)
    assign(paste("z1_",case,sep=""),z1)
}


#======#
# Plot #
#======#

# Plot params
beta_lim   = c(0,10)
z_lim      = c(0,30)  # km
zvec       = 2:(nz-1)
cex        = 1.5

pdf("~/Dropbox/17rad_cooling2/plots/beta.pdf",width=10,height=10,bg="white")
par(mfrow=c(2,2), mar=c(5,5,5,7))

for (k_case in 1:length(cases)){
	case	= cases[k_case]
	panel   = panels[k_case]
	casename= casenames[k_case]
    beta    = eval(as.name(paste("beta_",case,sep="")))
    beta1   = eval(as.name(paste("beta1_",case,sep="")))
    k       = eval(as.name(paste("k_",case,sep="")))
    nk	    = length(k) 
    z1      = eval(as.name(paste("z1_",case,sep="")))
    my.image.plot(1e-2*k,1e-3*z[zvec],beta,
	     ylim = z_lim,
	     zlim = beta_lim,
         xlab = "k (cm^-1)",
         ylab = "z (km)",
         main = bquote("("*.(panel)*")"~~frac(d~ln~tau[k],d~ln~p)*", "~.(casename)),
         cex.lab  = cex,
         cex.axis = cex,
         cex.main = cex,
         cex.legend = cex
         )
     ncoarse   = 10  # every 10 cm^-1
     z1_coarse = coarse_grain(z1,ncoarse)  
     k_coarse  = coarse_grain(k,ncoarse)
     points(1e-2*k_coarse,1e-3*z1_coarse,type="l",lty="dashed")

	 # beta_av
	 beta1_av =  round(mean(beta1,na.rm=TRUE),digits=1)
	 beta1_sd =  round(sd(beta1,na.rm=TRUE),digits=1)
	 text(0.8*max(1e-2*k),0.9*z_lim[2],bquote(beta[av]==.(beta1_av)%+-%.(beta1_sd)),cex=cex,col="white")
}
dev.off()

