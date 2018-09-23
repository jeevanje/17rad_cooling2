Rtoolsdir = "~/Dropbox/Rtools/"

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))

# params
cases	   = c("h2o_only_no_cont","co2_only_simple_atm")
Ncases     = length(cases)
k_lims	   = list(c(100,1500),c(500,850)) # cm^-1 
k_lab	   = expression(k~~"("*c*m^{-1}*")")
kappa_lab  = expression(kappa[ref]~~"("*m^2/kg*")")
kappa_lim  = c(1e-5,1e3) 
cex	       = 1.5
cex_leg    = 2
lwd        = 2
ncoarse    = 10
pref	   = ps
legend     = c(expression(dk==1~cm^{-1}),
			   expression(10~cm^{-1}~~"bins"),
			   "linear fit")
ltyvec     = c("dotted","dashed","solid")
colvec     = c("gray","darkgray","black")
legendposs = c("bottomleft","bottom")

# Band approximation
band_data  = list()
kj_gases   = list(c(1.5e4,1e5,1.45e5),c(500e2,667.5e2,850e2)) 

# PDF
file = "~/Dropbox/17rad_cooling2/plots/kappa.pdf"
pdf(file,width=10,height=5,bg="white")
par(mfrow=c(1,2),mar=c(5,5,5,3))

for (n in 1:Ncases){
	case    = cases[n]
    gas     = substr(case,1,3)	
	m_mol   = eval(as.name(paste("m_",gas,sep="")))
	k_lim   = k_lims[[n]]
	kj_vals = kj_gases[[n]]
	Nj      = length(kj_vals)

	#=======#
	# Data  #
	#=======#
	
	ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	p       = ncvar_get(nc,"p")
	dk      = k[2]-k[1]
	nk      = length(k)
	sigma   = ncvar_get(nc,"sigma") # m^2/molec.
	mref    = which.min(abs(p-ps))
	kappa   = as.array(sigma[ ,mref])*N_avo/m_mol     # m^2/kg, Tref=Ts
	
	#=========#
	# Process #
	#=========#
	
	kappa_raw    = kappa
	kappa_coarse = exp(coarse_grain(log(kappa),ncoarse))
	k_coarse 	 = coarse_grain(k,ncoarse)	
	ij_vals		 = numeric(Nj)
	for (j in 1:Nj){
		ij_vals[j] = which.min(abs(k_coarse-kj_vals[j]))
	}
	band_data_gas= list()
	for (j in 1:2){		
		ivec      = ij_vals[j]:ij_vals[j+1]
		k_fit     = k_coarse[ivec]
		fit		  = lm(log(kappa_coarse[ivec]) ~ k_fit)	
		lk     	  = -(fit$coefficients[[2]])^-1
		kappa_fit = exp(fit$fitted.values)
		
		data	  = list()
		for (var in c("k_fit","kappa_fit","lk","ivec")){
			data[[var]] <- abs(eval(as.name(var)))
		}	
		band_data_gas[[j]] <- data
	}
	band_data[[gas]] <- band_data_gas
	
	#=======#
	# Plot  #
	#=======#

	plot(1,type="n",
		xlim 	 = k_lim, 
		ylim	 = kappa_lim,
		xlab 	 = k_lab,
		ylab     = kappa_lab,
		log      = "y",
		main	 = paste(gas," absorption spectrum",sep=""),
		cex.lab  = cex,
		cex.axis = cex,
		cex.main = cex)
	# 1 cm^-1 plot
	points(1e-2*k,kappa_raw,type="l",lwd=1,lty=ltyvec[1],col=colvec[1])	
	# 10 cm^-1 plot
	points(1e-2*k_coarse,kappa_coarse,type="l",
		  lwd=4,lty=ltyvec[2],col=colvec[2])	

	# plot fits 
	for (j in 1:2){
		data = band_data_gas[[j]]	
		for (var in c("k_fit","kappa_fit","lk")){
			assign(var,data[[var]])
		}				
		points(1e-2*k_fit,kappa_fit,type="l",lwd=2,
		       lty=ltyvec[3],col=colvec[3])	
	} # bands
	legend(legendposs[n],legend=legend,lty=ltyvec,col=colvec,lwd=2.5,cex=1)
} # case
dev.off()
save(band_data,ncoarse,pref,k_coarse,file="~/Dropbox/17rad_cooling2/data/band_data.Rdata")