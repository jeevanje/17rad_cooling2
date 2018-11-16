Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")
library(ncdf4)
library(fields)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source("src/plot_params.R")
source("src/cts_functions.R")


#=======#
# Data  #
#=======#

fields    = c("coo2d","coo1d")
cases     = c("h2o_only_no_cont","h2o_only_simple_atm","h2o_co2_simple")
casenames = c("base","w/cont.","w/cont. + co2")

for (case in cases){
	ncpath  = paste("data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk      = k[2]-k[1]
	p       = ncvar_get(nc,"p")
	coo2d   = ncvar_get(nc,"coo")/abs(coo2d_fac)         # convert to SI 
	
	np      = length(p) 
	p_s     = rfm_lni2s(p)
	np_s    = length(p_s)
	coo2d_s = rfm_i2s(coo2d)
	coo1d_s = apply(coo2d_s,2,sum)*dk  # W/m^2/Pa
	k_coarse= coarse_grain(k,ncoarse)
	for (field in fields){
		assign(paste(field,"_",case,sep=""),eval(as.name(paste(field,"_s",sep=""))))
	}
}

#=======#
# Plot  #
#=======#

coo1d_lim = c(-6,0)   # K/day
cex = 2.5
lwd = 2.5
ltyvec = c("dotted","dashed","solid")
file = "plots/cont_co2_effects.pdf"
panelvec = c("(a)","(a)","(b)")
pdf(file,width=17,height=5,bg="white")
par(mfrow=c(1,3),mar=c(5,6,5,12))
for (k_case in 2:3){
	case  = cases[k_case]
	main  = bquote(.(panelvec[k_case])~~H[k]~~.(casenames[k_case]))
	coo2d = eval(as.name(paste("coo2d_",case,sep="")))
	coo2d_coarse = coarse_grain(coo2d,ncoarse)
	plot_coo2d(k_coarse,p_s,coo2d_coarse,main)
}
plot(1,type="n", 
     xlim     = coo1d_lim, 
     ylim     = p_lim,
     xlab     = "H (K/day)",    
     ylab     = plab,
     main     = "(c)  Heating rate",
     cex.lab  = cex,
     cex.axis = cex,
     cex.main = cex)
for (k_case in 1:length(cases)){
	case  = cases[k_case]
	lty	  = ltyvec[k_case]
	coo1d = eval(as.name(paste("coo1d_",case,sep="")))
	points(coo1d_fac*coo1d[-np_s],1e-2*p_s[-np_s],type="l",lwd=lwd,lty=lty,cex=cex)
}
legend(-6.2,50,legend=casenames,lty=ltyvec,lwd=lwd,cex=2)
dev.off()

