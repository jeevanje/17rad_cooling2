Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source("./src_paper/plot_params.R")
source("./src_paper/cts_functions.R")
load("./data_paper/band_params.Rdata")  # includes ncoarse

# Data params
case = "h2o_noctm_Ts300_rh0.75_gamma7"

#=======#
# Data  #
#=======#
gas     = "h2o"
ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk      = k[2]-k[1]
p       = ncvar_get(nc,"p")
np      = length(p)
z       = ncvar_get(nc,"z")
tabs    = ncvar_get(nc,"tabs")
Ts	    = tabs[1]
Tstrat  = tabs[np]
Tav     = 0.5*(Ts+Tstrat)
q_h2o   = ncvar_get(nc,"q_h2o")
q_co2   = ncvar_get(nc,"q_co2")
Gamma   = -(diff(tabs)/diff(z))[1]
RH      = (q_h2o/qsat(tabs,p))[1]	
nk      = length(k)
coo2d   = ncvar_get(nc,"coo")/abs(coo2d_fac)  # convert to SI 
opt     = D*ncvar_get(nc,"opt")
flx     = ncvar_get(nc,"flx")*1e-2			  # pppf units

# coo, cts
WVP        = calc_WVP0(tabs,RH,Gamma)*exp(-L/Rv/tabs)  # i lev
tau_ssm = kappa_simple_h2o%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*WVP)
p_s      = rfm_lni2s(p)
z_s      = rfm_i2s(z)
tabs_s   = rfm_i2s(tabs)
coo2d_s  = rfm_i2s(coo2d)
#coo2d_s  = (flx[ ,-np]-flx[ ,-1])/(rep(1,times=nk)%o%diff(p))
cts2d_s  = calc_cts2d(k,p,tabs_s,opt)
k_fit    = eval(as.name(paste("k_fit_",gas,sep="")))
dk_fit   = eval(as.name(paste("dk_fit_",gas,sep="")))
cts_ssm2d = calc_cts2d(k_fit,p,tabs_s,tau_ssm)
coo1d    = apply(coo2d_s,2,sum)*dk
coo1d_k150 = apply(coo2d_s[k>150e2,],2,sum)*dk
cts1d    = apply(cts2d_s,2,sum)*dk
cts_ssm1d = apply(cts_ssm2d,2,sum)*dk_fit	
cts_simple_bands = calc_cts_simple_bands(gas,p_s,tabs_s,Ts,Gamma,RH)	
cts_simple1d =  apply(cts_simple_bands,1,sum,na.rm=TRUE)
cts_simple1d[is.na(cts_simple1d)] <- 0


#==============#
# Plot prelims #
#==============#

# params
coo1d_lim  = c(-3,0)   # K/day
legendpos  = "bottomright"
cex	       = 1.25
cex_leg    = 1
lwd        = 2.5
fields1d   = c("coo1d","cts1d","cts_ssm1d","cts_simple1d","coo1d_k150")
ltyvec     = c("solid","solid","dashed","solid","dashed")
colvec     = c("black","gray","red","red","black")
legendvec  = c("RFM","RFM_CTS","SSM2D","SSM1D","RFM_k150")
fieldvec   = c(1,3,4)

# y coord
p_lim      = c(1000,0) # hPa
p_lab      = "p (hPa)"
p_fac	   = 1e-2        # SI to hPa
ycoord     = "p"
for (var in c("lim","lab","fac")){
	assign(paste("y_",var,sep=''),eval(as.name(paste(ycoord,"_",var,sep=''))))
}

plot_H_frame = function(main,coo1d_lim){
	plot(1,type="n", 
	   xlim     = coo1d_lim, 
	   ylim     = y_lim,
	   xlab     = Hlab,
	   ylab     = y_lab,
	   main     = main,
	   cex.lab  = cex,
	   cex.axis = cex,
	   cex.main = cex)	
}

#=======#
# Plot  #
#=======#
	
# PDF
file = "~/Dropbox/17rad_cooling2/plots_paper/H1d_validation.pdf"
pdf(file,width=5,height=5,bg="white")
par(mfrow=c(1,1),mar=c(5,5,5,3))

y_s = eval(as.name(paste(ycoord,"s",sep="_")))
main = "Heating rates for BASE"
plot_H_frame(main,coo1d_lim)
for (n_field in fieldvec){
    lty   = ltyvec[n_field]
    col   = colvec[n_field]
    field = eval(as.name(fields1d[n_field]))
    points(coo1d_fac*field,y_fac*y_s,type="l",lwd=lwd,lty=lty,col=col)  
}
legend(legendpos,legendvec[fieldvec],lty=ltyvec[fieldvec],
	  lwd=2,cex=cex_leg,col=colvec[fieldvec])
dev.off()	

