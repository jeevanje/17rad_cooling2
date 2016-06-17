library(ncdf)
library(fields)
source("~/Dropbox/Rtools/gray_model.R")

sst = 300
Ttp = 200
qvpert_list = c(0,-0.2,-0.4,-0.6,-0.8)
N = length(qvpert_list)
Qnet_vec = numeric(N)
taus_vec = numeric(N)

#===============#
# Get  data     #
#===============#

for (i in 1:N) {
    pert = qvpert_list[i]
    datadir=paste("~/Dropbox/rad_cooling2/data/rce",sst,"_h2o_only_nofusion_qvpert_",pert,sep="")
    ncpath<-paste(datadir,"/data/verticalstats.nc",sep="")
    open.ncdf(ncpath)->nc
    if (pert == 0) {
        get.var.ncdf(nc,"z")-> z    
	get.var.ncdf(nc,"tabs")->tabs  
	nz  = length(z)  
	ktp = which.min(abs(tabs-Ttp))
	}
    get.var.ncdf(nc,"qv")-> qv    
    get.var.ncdf(nc,"rho")-> rho    
    get.var.ncdf(nc,"lwup")-> lwup    
    get.var.ncdf(nc,"lwdown")-> lwdown
    get.var.ncdf(nc,"swup")-> swup    
    get.var.ncdf(nc,"swdown")-> swdown
    Fnet = lwup + swup - lwdown - swdown
    lapse = -partialder_i2s(3,z,s2i(3,z,tabs))
    pptfnet  = -partialder_i2s(3,z,Fnet)/lapse
    Qnet_vec[i] = Fnet[nz]-Fnet[1]
    assign(paste("Fnet_",pert,sep=""),Fnet)        
    assign(paste("pptfnet_",pert,sep=""),pptfnet)    
    assign(paste("rhov_",pert,sep=""),rho*qv)    
    close.ncdf(nc)
    }


#========#
# Plots  #
#========#

zvec=1:ktp
pptflim=c(0,5)
tabslim=c(sst,200)
cex=1.5
colvec = rev(tim.colors(N))
gray_col = "dimgrey"
rrtm_col = "red"
olr_col  = "black"
Qnet_col  = "black"
 
pdf(file=paste("../plots/pptfnet_varqv_300.pdf",sep=""),width=10,height=5,bg="white")
par(mfrow=c(1,2),mar=c(5,5,5,3))

# Panel 1 -- RRTM Qnet 
xlab = expression(q[v]~~"rescaling factor")
olrlim = range(Qnet_vec)
plot(1,type="n",
	ylab=expression(Q[net]~~"("~W/m^2~")"),
	xlab = xlab,
	main = "(a) Net cooling, RRTM",
	xlim = rev(range(1+qvpert_list)),
	ylim = olrlim,
	lty = "dashed",
	cex = cex,
	cex.main=cex,
	cex.lab=cex,
	cex.axis=cex)
points(1+qvpert_list,Qnet_vec,type="b",lty="dashed",
	pch=16,col=Qnet_col,cex=cex)

# Panel 2 -- RRTM pptfnet profiles
plot(1,type="n",xlim=pptflim,ylim=tabslim,
	xlab=expression(-partialdiff[T]*F[net]~~"("~W/m^2/K~")"),
	ylab="Temperature (K)",
	main = "(b) Net flux divergence, RRTM",
	cex.main=cex,
	cex.lab=cex,
	cex.axis=cex)

for (i in 1:N){
    pert =  qvpert_list[i]
    pptfnet = eval(as.name(paste("pptfnet_",pert,sep="")))
    points(-pptfnet[zvec],tabs[zvec],type="l",lwd=2,col=colvec[i])
    }

legend("topright",paste((1+qvpert_list),"qv",sep=" "),
	lty="solid",col=colvec,cex=1.25 )


dev.off()
