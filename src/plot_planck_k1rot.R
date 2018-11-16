Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("./src/plot_params.R")
source("./src/cts_functions.R")

# params
case      = "h2o_only_no_cont"
RH        = 0.75
Gamma     = 7e-3            # K/km

# data
ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k") 
p       = ncvar_get(nc,"p") 
tabs    = ncvar_get(nc,"tabs")
Ts      = tabs[1]
np      = length(p) 
pvec	= np:1
k1_rot  = calc_k1("h2o",p,tabs,Ts,Gamma,RH)[[1]]
beta    = 1 + L/Rv/tabs*Gamma*Rd/g

# plot params
cex     = 1.25
cex_text = 1.2
lwd 	= 2
pmin    = p[min(which(is.na(k1_rot)))]
plim	= 1e-2*rev(range(p))

pdf("./plots/planck_k1rot.pdf",width=10,height=5)
layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1,1))

# Planck
par(mar=c(5,5,4,2))
my.image.plot(1e-2*k,1e-2*p[pvec],1e2*pi*outer(k,tabs[pvec],planck_k),
	xlim     = 1e-2*range(k),
	ylim	 = plim,
	xlab     = klab,
	ylab     = plab,
    main     = expression("Planck function"~~ "("*W/m^2/cm^{-1}*")"),
	cex.main = cex,
	cex.axis = cex,
	cex.lab  = cex,
	cex.legend = 1.25)
points(1e-2*k1_rot,1e-2*p,type="l",lty="dashed",lwd=lwd+0.5)
ptext = 3e4 # hPa
ktext = which.min(abs(p-ptext))
text(500,1e-2*ptext,expression(k[1*rot]),cex=1.5)

# beta/p
par(mar=c(5,5.5,4,1.5))
plot(beta[p>pmin]/beta[1]*ps/p[p>pmin],1e-2*p[p>pmin],
	type="l",xlab="",ylab=plab,lwd=lwd,
	xlim	 = c(2e-1,5),
	ylim     = plim,
	log      = "x",
	main     = "Planck emission vs. emissivity gradient",
	cex.axis = cex,
	cex.main = cex,
	cex.lab  = cex)
points(planck_k(k1_rot,tabs)/planck_k(k1_rot[1],tabs[1]),1e-2*p,
	type="l",lwd=lwd,lty="solid")
abline(v=1,lty="dotted",lwd=2)  
text(0.5,180,expression(B^{"*"}*"("*k[1*rot]*","*T*")"),cex=cex_text)	
text(4,200,expression(beta^{"*"}/p^{"*"}),cex=cex_text)
dev.off()
