Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source("./src/plot_params.R")
source("./src/cts_functions.R")
load("./data/band_params.Rdata")  # includes ncoarse

# Data params
case = "h2o_only_no_cont"

#=======#
# Data  #
#=======#

ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
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
q       = ncvar_get(nc,"q_h2o")
Gamma   = -(diff(tabs)/diff(z))[1]
RH      = (q/qsat(tabs,p))[1]	
nk      = length(k)
coo2d   = ncvar_get(nc,"coo")/abs(coo2d_fac)  # convert to SI 
opt     = D*ncvar_get(nc,"opt")

# coo, cts
WVP      = calc_WVP0(tabs,RH,Gamma)*exp(-L/Rv/tabs)  # i lev
tau_simple = kappa_simple%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*WVP)
p_s      = rfm_lni2s(p)
z_s      = rfm_i2s(z)
tabs_s   = rfm_i2s(tabs)
coo2d_s  = rfm_i2s(coo2d)
cts2d_s  = calc_cts2d(k,p,tabs_s,opt)
cts_simple2d = calc_cts2d(k_fit,p,tabs_s,tau_simple)
coo1d    = apply(coo2d_s,2,sum)*dk
cts1d    = apply(cts2d_s,2,sum)*dk
cts_simple1d = apply(cts_simple2d,2,sum)*dk_fit	
cts_verysimple1d = calc_cts_verysimple("h2o",p_s,tabs_s,Ts,Gamma,RH)	


#==============#
# Plot prelims #
#==============#

# params
coo1d_lim  = c(-3,0)   # K/day
cex	       = 1.25
cex_leg    = 1
lwd        = 2.5
fields1d   = c("coo1d","cts1d","cts_verysimple1d")
ltyvec     = c("solid","solid","dashed")
colvec     = c("black","gray","black")
legendvec  = c(expression(H[RFM]),expression(C*T*S[RFM]),"Theory")

# y coord
p_lim      = c(1000,0) # hPa
p_lab      = "p (hPa)"
p_fac	   = 1e-2        # SI to hPa
ycoord     = "p"
for (var in c("lim","lab","fac")){
	assign(paste("y_",var,sep=''),eval(as.name(paste(ycoord,"_",var,sep=''))))
}

plot_H_frame = function(main){
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
file = "~/Dropbox/17rad_cooling2/plots/H1d_rfm_cts_theory.pdf"
pdf(file,width=5,height=5,bg="white")
par(mar=c(5,5,5,3))

y_s        = eval(as.name(paste(ycoord,"s",sep="_")))
plot_H_frame("Heating rates, RFM + theory")
for (n in 1:length(fields1d)){
    lty   = ltyvec[n]
    col   = colvec[n]
    field = eval(as.name(fields1d[n]))
    points(coo1d_fac*field,y_fac*y_s,type="l",lwd=lwd,lty=lty,col=col)  
}
#points(coo1d_fac*cts_simple1d,y_fac*y_s,type="l",lwd=lwd-1,lty="dotted",col="black")  
legend("bottomright",legendvec,lty=ltyvec,lwd=2,cex=cex_leg,col=colvec)
dev.off()	

