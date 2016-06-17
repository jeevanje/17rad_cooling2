library(ncdf)
library(fields)
source("~/Dropbox/Rtools/gray_model.R")

sst = 300
Ttp = 200
qvpert_list = c(0,-0.2,-0.4,-0.6,-0.8)
N = length(qvpert_list)
Qlw_vec = numeric(N)
OLR_vec = numeric(N)
Qlw_gray_vec = numeric(N)
OLR_gray_vec = numeric(N)
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
    Flw = lwup - lwdown
    lapse = -partialder_i2s(3,z,s2i(3,z,tabs))
    pptflw  = -partialder_i2s(3,z,Flw)/lapse
    Qlw_vec[i] = Flw[nz]-Flw[1]
    OLR_vec[i] = Flw[nz]
    assign(paste("Flw_",pert,sep=""),Flw)        
    assign(paste("pptflw_",pert,sep=""),pptflw)    
    assign(paste("rhov_",pert,sep=""),rho*qv)    
    close.ncdf(nc)
    }


#========#
# Plots  #
#========#

zvec=1:ktp
pptflim=c(0,12)
tabslim=c(sst,200)
cex=1.5
colvec = rev(tim.colors(N))
kappa = tune_kappa(z,rhov_0,tabs,sst,OLR_vec[1])  # m2/kg
gray_col = "dimgrey"
rrtm_col = "red"
olr_col  = "black"
Qlw_col  = "black"
olrlim = c(150,400)
 
pdf(file=paste("../plots/olr_qlw_window.pdf",sep=""),width=10,height=10,bg="white")
par(mfrow=c(2,2),mar=c(5,5,5,3))

# Panel 1 -- RRTM Qlw and OLR
xlab = expression(q[v]~~"rescaling factor")
#olrlim = range(c(OLR_vec,Qlw_vec))
plot(1,type="n",
	ylab=expression(Q[LW]*", OLR"~~"("~W/m^2~")"),
	xlab = xlab,
	main = "(a) LW cooling and OLR, RRTM",
	xlim = rev(range(1+qvpert_list)),
	ylim = olrlim,
	lty = "dashed",
	cex = cex,
	cex.main=cex,
	cex.lab=cex,
	cex.axis=cex)
points(1+qvpert_list,Qlw_vec,type="b",lty="dashed",
	pch=16,col=Qlw_col,cex=cex)
points(1+qvpert_list,OLR_vec,type="b",lty="dashed",
	pch=16,col=olr_col,cex=cex)
text(0.4,225,expression(Q[LW]),cex=cex)
text(0.4,325,expression("OLR"),cex=cex)

# Panel 2 -- RRTM pptf profiles
plot(1,type="n",xlim=pptflim,ylim=tabslim,
	xlab=expression(-partialdiff[T]*F~~"("~W/m^2/K~")"),
	ylab="Temperature (K)",
	main = "(b) LW flux divergence, RRTM",
	cex.main=cex,
	cex.lab=cex,
	cex.axis=cex)

for (i in 1:N){
    pert =  qvpert_list[i]
    pptflw = eval(as.name(paste("pptflw_",pert,sep="")))
    points(-pptflw[zvec],tabs[zvec],type="l",lwd=2,col=colvec[i])
    }

legend("topright",paste((1+qvpert_list),"qv",sep=" "),
	lty="solid",col=colvec,cex=1.25 )

# Panel 3 -- Gray Qlw and OLR
for (i in 1:N){
	pert =  qvpert_list[i]
    rhov = eval(as.name(paste("rhov_",pert,sep="")))
	wvp  = compute_wvp(z,rhov)
	taus = kappa*wvp[1]
    U    = compute_Ugray(kappa,z,rhov,tabs,sst)
    D    = compute_Dgray(kappa,z,rhov,tabs)
    F    = U-D
    pptflw = compute_pptf(kappa,z,U,D,tabs,rhov, lapse)
 	taus_vec[i]	    = taus   
    Qlw_gray_vec[i] = F[nz]-F[1]
    OLR_gray_vec[i] = F[nz]
    }
#olrlim = range(c(OLR_gray_vec,Qlw_gray_vec))
xlab = expression(q[v]~~"rescaling factor")
plot(1,type="n",
    ylab=expression(Q[LW]*", OLR"~~"("~W/m^2~")"),
    xlab = xlab,
    main = "LW cooling and OLR, gray",
    xlim = rev(range(1+qvpert_list)),
    ylim = olrlim,
    lty = "dashed",
    cex = cex,
    cex.main=cex,
    cex.lab=cex,
    cex.axis=cex)
points(1+qvpert_list,Qlw_gray_vec,type="b",lty="dashed",
    pch=16,col=Qlw_col,cex=cex)
points(1+qvpert_list,OLR_gray_vec,type="b",lty="dashed",
    pch=16,col=olr_col,cex=cex)
text(0.4,225,expression(Q[LW]),cex=cex)
text(0.4,325,expression("OLR"),cex=cex)

# Panel 2 --  pptf profiles
plot(1,type="n",xlim=pptflim,ylim=tabslim,
    xlab=expression(-partialdiff[T]*F~~"("~W/m^2/K~")"),
    ylab="Temperature (K)",
    main = "LW flux divergence, gray",
    cex.main=cex,
    cex.lab=cex,
    cex.axis=cex)   
for (i in 1:N){
    pert =  qvpert_list[i]
    rhov = eval(as.name(paste("rhov_",pert,sep="")))
    U    = compute_Ugray(kappa,z,rhov,tabs,sst)
    D    = compute_Dgray(kappa,z,rhov,tabs)
    pptflw = compute_pptf(kappa,z,U,D,tabs,rhov, lapse)
    points(-pptflw[zvec],tabs[zvec],type="l",lwd=2,col=colvec[i])
    }

legend("topright",paste((1+qvpert_list),"qv, ",expression(tau_s),"=",round(taus_vec,digits=1),sep=" "),
                lty="solid",col=colvec,cex=1.25 )
                


dev.off()
