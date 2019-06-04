Rtoolsdir = "~/Dropbox/Rtools/"

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))

# params
cases      = c("h2o_noctm_Ts300_rh0.75_gamma7","co2_noctm_Ts300_rh0.75_gamma7")
Ncases     = length(cases)
pref	   = ps

# Band approximation
band_data  = list()
kj_gases   = list(c(1.5e4,1e5,1.45e5),c(500e2,667.5e2,850e2)) 

for (n in 1:Ncases){
	case    = cases[n]
    gas     = substr(case,1,3)	
	m_mol   = eval(as.name(paste("m_",gas,sep="")))
	kj_vals = kj_gases[[n]]
	Nj      = length(kj_vals)
	sigmavar = paste("sigma",gas,sep="_")

	#=======#
	# Data  #
	#=======#
	
	ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	p       = ncvar_get(nc,"p")
	dk      = k[2]-k[1]
	nk      = length(k)
	sigma   = ncvar_get(nc,sigmavar) # m^2/molec.
	mref    = which.min(abs(p-ps))
	kappa   = as.array(sigma[ ,mref])*N_avo/m_mol     # m^2/kg, Tref=Ts
	ncoarse = 1000/dk
	
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
} # case
save(band_data,ncoarse,pref,k_coarse,file="~/Dropbox/17rad_cooling2/data_paper/band_data.Rdata")