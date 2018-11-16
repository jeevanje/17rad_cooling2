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
Tslist = seq(270,300,by=10)
cases  = paste("h2o_only_no_cont_",Tslist,"K",sep="")
Ncase = length(cases)

#=======#
# Data  #
#=======#

for (case in cases){
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
	dtabsdp = diff(tabs)/diff(p) 

	# coo, cts
	p_s      = rfm_lni2s(p)
	z_s      = rfm_i2s(z)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	coo1d    = apply(coo2d_s,2,sum)*dk
	cts_verysimple1d = calc_cts_verysimple("h2o",p_s,tabs_s,Ts,Gamma,RH)	

	vars     = c("coo1d","cts_verysimple1d","z_s",
				 "p_s","tabs_s","Ts","Gamma","RH","dtabsdp")
	for (var in vars){
		assign(paste(var,Ts,sep="_"),eval(as.name(var)))	
	}
}


#==============#
# Plot prelims #
#==============#

# common params
coo1d_lim  = c(-5,0)   # W/m^2/K
tabs_lim   = c(300,200) # K
pptflab    = expression(partialdiff[T]*F~~"("*W/m^2/K*")")
tabs_lab   = "Temperature (K)"
cex	       = 1.5
cex_leg    = 1.25
lwd        = 2.5
colvec     = tim.colors(Ncase)

plot_pptf_frame = function(main){
	plot(1,type="n", 
	   xlim     = coo1d_lim, 
	   ylim     = tabs_lim,
	   xlab     = pptflab,
	   ylab     = tabs_lab,
	   main     = main,
	   cex.lab  = cex,
	   cex.axis = cex,
	   cex.main = cex)	
}

#=======#
# Plot  #
#=======#
	
# PDF
file = "~/Dropbox/17rad_cooling2/plots/pptf_tinv.pdf"
pdf(file,width=9,height=5,bg="white")
par(mfrow=c(1,2),mar=c(5,5,4,3))


# RFM
plot_pptf_frame("RFM")
for (n in 1:Ncase){
	Ts    = Tslist[n] 
    col   = colvec[n]
    tabs_s= eval(as.name(paste("tabs_s",Ts,sep="_")))
    coo1d = eval(as.name(paste("coo1d",Ts,sep="_")))
    points(-coo1d/dtabsdp,tabs_s,type="l",lwd=lwd,lty="solid",col=col)  
}
legend("topleft",legend=Tslist,lwd=2,cex=cex_leg,col=colvec,
	title=expression(T[s]~~"(K)"))

# Theory
plot_pptf_frame("Theory")
for (n in 1:Ncase){
	Ts    = Tslist[n] 
    col   = colvec[n]
    tabs_s= eval(as.name(paste("tabs_s",Ts,sep="_")))
    cts_verysimple1d = eval(as.name(paste("cts_verysimple1d",Ts,sep="_")))
    points(-cts_verysimple1d/dtabsdp,tabs_s,type="l",lwd=lwd,
    		lty="solid",col=col)  
}
legend("topleft",legend=Tslist,lwd=2,cex=cex_leg,col=colvec,
	title=expression(T[s]~~"(K)"))

dev.off()	

