Rtoolsdir = "~/Dropbox/Rtools/"

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"my_image_plot.R",sep=""))
source("plot_params.R")

# Data params
cases	  = c("h2o_only_no_cont","co2_only_simple_atm")
#cases     = c("h2o_only_no_cont_280K")
Ncase     = length(cases)
D         = 1.5  # diffusion parameter
Nsec      = 3600*24   # sec/day
coo_fac   = -g/Cp*Nsec*1e2  # pppf to K/day/cm^-1 heating
cts_fac   = coo_fac         # pppf to K/day/cm^-1 heating
lntau_fac = 1 
fac1d	  = -g/Cp*Nsec      # pppf to K/day

# Plot params
# y coord
k_lim	   = c(100,1500) # cm^-1  
k_lab	   = expression(k~~"("*c*m^{-1}*")")
p_lim      = c(1000,1) # hPa
p_lab      = "p (hPa)"
p_fac	   = 1e-2        # SI to hPa
z_lim      = c(1,30)     # km
z_lab      = "z (km)"
z_fac	   = 1e-3	     # SI to km
ycoord     = "p"
for (var in c("lim","lab","fac")){
	assign(paste("y_",var,sep=''),eval(as.name(paste(ycoord,"_",var,sep=''))))
}

# fields
fields 	    = c("coo","cts","lntau")
fields1d    = c("coo","cts")
coo_lim     = c(-0.010,0.010) # K/day/cm^-1
cts_lim     = coo_lim
lntau_lim   = c(-20,20)
coo_units   = "(K/day/cm^-1)"
cts_units   = coo_units
lntau_units = ""
coo1d_lims  = list(c(-2.5,0),c(-0.75,0.1))   # K/day

kappa_lab   = expression(kappa[ref]~~"("*m^2/kg*")")
kappa_lim   = c(1e-4,1e4) 
Ncol	    = 2^6
colvec      = c(two.colors(round(-Ncol*coo_lim[1]/(coo_lim[2]-coo_lim[1])),
			   			  start='darkblue',middle='blue',end='white'), 
               two.colors(round(Ncol*coo_lim[2]/(coo_lim[2]-coo_lim[1])),
               			  start='white',middle='red',end='darkred') )
cex	       = 2.25
cex_leg    = 2
lwd        = 2
ltyvec     = c("solid","dashed")

nrun   	   = 10  # every 10 cm^-1
ncoarse	   = 10

for (i in 1:Ncase){
	case  = cases[i]

	# PDF
	file = paste("~/Dropbox/17rad_cooling2/plots/coo_tau_kappa_",case,".pdf",sep="")
	pdf(file,width=16,height=9,bg="white")
	par(mfrow=c(2,3),mar=c(5,7,5,9))

	#=======#
	# Data  #
	#=======#

	gas     = substr(case,1,3)	
	gasname = gasnames[i]
	m_mol   = eval(as.name(paste("m_",gas,sep="")))	
	ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk      = k[2]-k[1]
	p       = ncvar_get(nc,"p")
	z       = ncvar_get(nc,"z")
	y  	    = eval(as.name(ycoord))
	tabs    = ncvar_get(nc,"tabs")
	q       = ncvar_get(nc,paste("q_",gas,sep=""))
	rho_gas = p/Rd/tabs*q
	path    = round(rfm_i2s(rho_gas)%*%diff(z),digits=1)
	np      = length(p)
	p_vec	= (np-1):1
	z_vec   = 1:(np-1)
	y_vec   = eval(as.name(paste(ycoord,"_vec",sep="")))
	nk      = length(k)
	coo2d   = ncvar_get(nc,"coo")/abs(coo_fac)         # convert to SI 
	opt     = D*ncvar_get(nc,"opt")
	flx2d   = ncvar_get(nc,"flx")*1e-2          # W/m^2/m^-1
	kappa   = ncvar_get(nc,"sigma",start=c(1,1),count=c(nk,1))*N_avo/m_mol # m^2/kg
	
	#=========#
	# Process #
	#=========#
	
	p_s      = rfm_lni2s(p)
	y_s      = rfm_i2s(y)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	opt_s    = rfm_lni2s(opt)
	lntau2d_s= log(opt_s)
	dtaudp   = (opt[ ,2:np]-opt[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	B_s      = outer(k,tabs_s,planck_k)
	
	cts2d_s  = pi*B_s*exp(-opt_s)*dtaudp  #s, p_s, W/m^2/Pa/m^-1, pppf units
	pppf2d   = (flx2d[ ,2:np]-flx2d[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	#coo2d_s  = -pppf2d   # W/m^2/Pa/m^-1
	coo1d    = apply(coo2d_s,2,sum)*dk  # W/m^2/Pa
	cts1d    = apply(cts2d_s,2,sum)*dk
	
    p1      = array(dim=c(nk))
    for (m in 1:nk){
        p1[m] = p_s[which.min(abs(opt_s[m,]-1))]
        #p1[m] = p_s[max(which(opt[m,]>=1))]
    }
	
	#=======#
	# Plot  #
	#=======#

	coo1d_lim   = coo1d_lims[[i]]

	# 1Dfields 
    plot(1,type="n", 
        xlim     = coo1d_lim, 
        ylim     = y_lim,
        xlab     = "H (K/day)",    
        ylab     = y_lab,
        main     = paste("(a) ",gasname," heating rate",sep=""),
        cex.lab  = cex,
        cex.axis = cex,
        cex.main = cex)
    for (n in 1:length(fields1d)){
        field = fields1d[n]
        lty   = ltyvec[n]
        var   = eval(as.name(paste(field,"1d",sep="")))
        points(fac1d*var,y_fac*y_s,type="l",lwd=2,lty=lty)  
     }
legend("bottomright",legend=c(expression(H),expression(H^CTS)),lwd=lwd,
			lty=c("solid","dashed"),cex=2)

	# H2O rot. band
	if (gas=="h2o"){
		i_rot   = 1:which.min(abs(k-1e5))
		coo_rot1d = apply(coo2d_s[i_rot,],2,sum)*dk  # W/m^2/Pa
		cts_rot1d = apply(cts2d_s[i_rot,],2,sum)*dk  # W/m^2/Pa
        #points(fac1d*cts_rot1d,y_fac*y_s,type="l",lwd=2,lty="dotted")  
	}

	# 2D fields
	for (k_field in 1:length(fields)){
		field	   = fields[k_field]
	    var        = eval(as.name(paste(field,"2d_s",sep="")))
   	    var_coarse = coarse_grain(var,ncoarse)
	    k_coarse   = coarse_grain(k,ncoarse)
	    zlim       = eval(as.name(paste(field,"_lim",sep="")))
	    zfac       = eval(as.name(paste(field,"_fac",sep="")))
	    units      = eval(as.name(paste(field,"_units",sep="")))
		if (field=="coo"){
			main = bquote("(b)"~~.(gasname)~~H[k]*"  ("~K/day/cm^{-1}~")")
		} else if (field=="cts"){
			main = bquote("(c)"~~.(gasname)~~H[k]^{CTS}*"  ("~K/day/cm^{-1}~")")
		} else if (field=="lntau"){
			main = bquote("(d)"~~.(gasname)~~ln~tau[k])
		}
my.image.plot(1e-2*k_coarse,y_fac*y_s[y_vec],zfac*var_coarse[ ,y_vec],
			#xlim = k_lim,
			ylim = y_lim,
			zlim = zlim,
			xlab = k_lab,
			ylab = y_lab,
			log  = "",
			main = main,
			#col  = colvec,
			cex.lab  = cex,
		    cex.axis = cex,
			cex.main = cex,
			cex.legend = cex_leg
		 )
	     p1_coarse = coarse_grain(p1,ncoarse)  
	     points(1e-2*k_coarse,1e-2*p1_coarse,type="l",lty="dashed",lwd=1)
	}  # 2D fields

	# Kappa
	ncoarse = ncoarse
	kappa_coarse = coarse_grain(kappa,ncoarse)
	plot(1e-2*k_coarse,kappa_coarse,type="l",lwd=lwd,
		#xlim 	 = k_lim, 
		ylim	 = kappa_lim,
		xlab 	 = k_lab,
		ylab     = kappa_lab,
		log      = "y",
		main	 = paste("(e) ",gasname," absorption spectrum",sep=""),
		cex.lab  = cex,
		cex.axis = cex,
		cex.main = cex)

	dev.off()	
}

