Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("src/plot_params.R")
load("data_paper/band_data.Rdata")
load("data_paper/band_params.Rdata")

# Data params
case_h2o  = "h2o_noctm_Ts300_rh0.75_gamma7"
case_co2  = "co2_noctm_Ts300_rh0.75_gamma7"
D         = 1.5  # diffusion parameter
RH		  = 0.75
Gamma	  = 7e-3
gases     = c("h2o","co2")
gas	      = gases[1]
case      = eval(as.name(paste("case",gas,sep="_")))

#=======#
# RFM   #
#=======#

gasname = gasnames[1]
m_mol   = eval(as.name(paste("m_",gas,sep="")))
sigmavar = paste("sigma",gas,sep="_")	
ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
nc      = nc_open(ncpath)
k       = ncvar_get(nc,"k")
dk	    = k[2]-k[1]
p       = ncvar_get(nc,"p")
p_s     = rfm_lni2s(p)
q       = ncvar_get(nc,"q_co2")
q_s     = rfm_i2s(q)
tabs    = ncvar_get(nc,"tabs")
tabs_s  = rfm_i2s(tabs)
Ts      = tabs[1]
Tav     = mean(range(tabs))
np      = length(p)
p_vec	= (np-1):1
nk      = length(k)
ncoarse	= 1000/dk
mref    = which.min(abs(p-pref))

k_rfm    = k
kappa_rfm= ncvar_get(nc,"sigma_h2o")*N_avo/m_h2o
tau_rfm  = D*ncvar_get(nc,"opt")
coo	    = ncvar_get(nc,"coo")/abs(coo2d_fac)         # convert to SI 
coo_rfm = rfm_i2s(coo)

#=========#
# Theory  #
#=========#

band_data_gas = band_data[[gas]]

for (j in 1:2){
   data = band_data_gas[[j]]       
   for (var in c("k_fit","kappa_fit","lk")){
       assign(paste(var,j,sep=""),data[[var]])
   }
}
k_fit     = c(k_fit1,k_fit2[-1])
dk_fit    = k_fit[2] - k_fit[1]
if (gas=="h2o"){
	WVP_s 		  = D*RH*esat(tabs_s)*Tav/L/Gamma  # Koll 2018, include D
	WVP 		  = D*RH*esat(tabs)*Tav/L/Gamma  
	kappa_theory  = 0.5*c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
	                   kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr)) #0.5 for temp scaling!                        
	tau_theory    = kappa_theory%o%(p/pref*WVP)    
	tau_theory_s  = kappa_theory%o%(p_s/pref*WVP_s)         
} else if (gas=="co2"){
	kappa_theory  = kappa_Q*exp(-abs(k_fit-k_Q)/l_Q)         
	tau_theory    = 0.5*kappa_theory%o%(q*p^2/(g*pref))
	tau_theory_s  = 0.5*kappa_theory%o%(q_s*p_s^2/(g*pref))
}
nk_fit  = length(k_fit)
dtaudp  = (tau_theory[ ,-1]-tau_theory[ ,-np])/(rep(1,times=nk_fit)%o%diff(p)) #s lev
B       = outer(k_fit,tabs_s,planck_k)
trans   = exp(-tau_theory_s)
coo_theory = pi*B*trans*dtaudp  # s levs
k_theory   = k_fit

#=======#
# Plot  #
#=======#

# Plot params
# y coord
p_lim      = c(1000,10) # hPa
p_fac	   = 1e-2        # SI to hPa
p_lab	   = plab
chk_lab    = expression(H[tilde(nu)]*"  ("*K/day/cm^{-1}*")")
coo_lim    = c(-0.0125,0) # K/day/cm^-1
lnopt_lim  = c(-15,5)
coo_units  = "(K/day/cm^-1)"
cex	       = 1.5
cex_leg    = 1
lwd        = 3
main	   = ""

methods    = c("rfm","theory","rfm_avg")
methodnames= c("RFM","SSM2D","RFM_avg")
cols       = c("black","red","black")
ltys       = c("solid","solid","dashed")
lwds       = c(2,2,4)
kvals      = c(5e4)  # m^-1
Qvals      = numeric(3)

# PDF
file = "~/Dropbox/17rad_cooling2/plots_paper/chk_vs_p.pdf"
pdf(file,width=5,height=6,bg="white")
par(mfrow=c(1,1),mar=c(5,6,3,3))
for (kval in kvals){
	i_kval_theory        = which.min(abs(k_theory-kval))
	kapparef_kval_theory = kappa_theory[i_kval_theory]
	i_kval_rfm1  	     = which.min(abs(k_rfm-kval))

	# kappavals
	kapparef_kval_rfm    = kappa_rfm[i_kval_rfm1,mref]
	i_kval_rfm = which.min(abs(kappa_rfm[(i_kval_rfm1-100):(i_kval_rfm1+100),mref]-kapparef_kval_theory)) + (i_kval_rfm1-100-1) 
	kapparef_kappaval = kappa_rfm[i_kval_rfm,mref]

	i_kvals_rfm 	     = (i_kval_rfm-50):(i_kval_rfm+50)
	
	# chk profiles
	plot(1,type="n",xlim=coo_lim,ylim=plim,
		 xlab = chk_lab,
		 ylab = plab,
		 main = main,
		 cex.lab  = cex,
		 cex.axis = cex,
		 cex.main = cex)
	for (n in 1:2){
		method    = methods[n]
		i_kval    = eval(as.name(paste("i_kval_",method,sep="")))
		col       = cols[n]
		lty       = ltys[n]
		lwd       = lwds[n]
		coo2dvals = eval(as.name(paste("coo_",method,sep="")))
		coovals   = coo2dvals[i_kval, ]
		Qvals[n]  = round(1e2*coovals%*%diff(p),digits=4) # W/m^2/cm-1
		points(coo2d_fac*coovals,1e-2*p_s,type="l",col=col,lty=lty,lwd=lwd)
	}
	coovals_avg = apply(coo_rfm[i_kvals_rfm,],2,mean)
	Qvals[3]    = round(1e2*coovals_avg%*%diff(p),digits=4) # W/m^2/cm-1
	points(coo2d_fac*coovals_avg,1e-2*p_s,type="l",col=cols[3],lty=ltys[3],lwd=lwds[3])
	legend("topleft",legend=methodnames,lty=ltys,col=cols,lwd=lwd[1],cex=1.1)

	#legend("bottomleft",legend=round(Qvals,digits=3),
	#		 title=expression(frac(C[p],g)*integral(H[tilde(nu)]*d*p,0,p[s])),
	#		 lty=ltys,col=cols,lwd=lwd[1],cex=1,bty="n")
	# Integrals agree to within 5%
	
	# # opt profiles
	# plot(1,type="n",xlim=lnopt_lim,ylim=plim,
		 # xlab = lntaulab,
		 # ylab = plab,
		 # main = main,
		 # cex.lab  = cex,
		 # cex.axis = cex,
		 # cex.main = cex)
	# for (k_method in 1:2){
		# method    = methods[k_method]
		# col       = cols[k_method]
		# lty       = ltys[k_method]
		# i_kval    = eval(as.name(paste("i_kval_",method,sep="")))
		# tau2dvals = eval(as.name(paste("tau_",method,sep="")))
		# tauvals   = tau2dvals[i_kval, ]
		# points(log(tauvals),1e-2*p,type="l",col=col,lty=lty,lwd=lwd)
	# }
} 
dev.off()	


