Rtoolsdir = "~/Dropbox/Rtools/"
#setwd("~/Dropbox/17rad_cooling2/src/")

library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"gray_model.R",sep=""))
source("plot_params.R")

# data params
cases = c("h2o_only_no_cont","co2_only_simple_atm")
D     = 1.5  # diffusion parameter
Nsec  = 3600*24   # sec/day
fac   = g/Cp*Nsec*1e2  # pppf, SI to K/day/cm^-1

# plot params
coo_lim = c(-0.005,0.015) # K/day/cm^-1
fields 	= c("coo","sx","ax","gex","cts","sum")
fieldnames = c(expression(H[k]),"SX","AX","GX","CTS","sum")
N_fields= length(fields)-1 # no sum
ltyvec	= c("solid",rep("dashed",time=5))
colvec  = c("black","forestgreen","orange","red","blue","black")
cex		= 2
plab    = "p (hPa)"
p_lim   = c(1000,0.1)

# PDF
file = "~/Dropbox/17rad_cooling2/plots/cooling_profiles.pdf"
pdf(file=file,width=12,height=10,bg="white")
par(mfrow=c(2,3),mar=c(5,5,5,4))

for (i in 1:2){
	case = cases[i]

	#=======#
	# Data  #
	#=======#
	
	# Get data
	ncpath  = paste("~/Dropbox/17rad_cooling2/data/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk		= k[2]-k[1]
	p       = ncvar_get(nc,"p")
	np      = length(p)
	tabs    = ncvar_get(nc,"tabs")
	nk		= length(k)
	coo2d   = ncvar_get(nc,"coo")/fac     # convert to SI 
	opt     = D*ncvar_get(nc,"opt")
	flx2d   = ncvar_get(nc,"flx")*1e-2    # W/m^2/m-1
	
	#=========#
	# Process #
	#=========#
	
	p_s	 	 = rfm_lni2s(p)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	opt_s	 = rfm_lni2s(opt)
	dtaudp2d = (opt[ ,2:np]-opt[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	B_i      = outer(k,tabs,planck_k)
	B_s      = outer(k,tabs_s,planck_k)
	opt_surf = opt[ ,1]%o%(rep(1,times=(np-1)))
	B_surf   = B_i[ ,1]%o%(rep(1,times=(np-1)))
	
	cts2d	  = pi*B_s*exp(-opt_s)*dtaudp2d  #s, z_s, W/m^2/Pa/m^-1
	gex2d     = -pi*(B_surf-B_s)*exp(-(opt_surf - opt_s))*dtaudp2d  #s, z_s
	pppf2d    = (flx2d[ ,2:np]-flx2d[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #slev
	
	#=======#
	# Plot  #
	#=======#
	gas     = substr(case,1,3)
	gasname = gasnames[i]

	p1vals  = 1e2*c(300,550,800)
	Np1     = length(p1vals)
	mvals	= numeric(Np1)
	kmin    = 3e4
	mmin    = min(which(k>kmin))

    p1vec   = numeric(nk)
    for (m in 1:nk){
        p1vec[m] = p[which.min(abs(opt[m,]-1))]
    }
	
	for (n in 1:Np1){
	    p1       = p1vals[n]
	    mvals[n] = which.min(abs(p1vec[(k>kmin)&(k<10e4)]-p1)) + mmin - 1
	}
	
	for (m in mvals){
	    kval    = k[m] 
	    tauvals = opt[m,1:(np-1)]  #i
	    Bvals   = B_s[m,] 	       #s
	    dtaudp  = dtaudp2d[m,]     #s
	    pppf    = pppf2d[m,]
	    coo   = -pppf 	# use simple difference rather than RFM's 3-stencil
	    cts   = cts2d[m,]
	    gex   = gex2d[m,]
	    taus  = round(tauvals[1],digits=1)
	
	    # gray_model.R
	    ex    = compute_ex(pi*Bvals,tauvals)  #B s levs, tau i levs
	    ax    = (ex$ax)*dtaudp
	    sx    = (ex$sx)*dtaudp
	    sum   = cts+gex+ax+sx
	
	    plot(1,type="n",
		    ylim = p_lim,
		    xlim = coo_lim,
		    xlab = "Cooling (K/day/cm^-1)",
		    ylab = plab,
		    main = bquote(.(gasname)*", k="~.(1e-2*kval)~c*m^{-1}*","~tau[s]==.(taus)),
		    cex.lab  = cex,
		    cex.axis = cex,
		    cex.main = cex
	    )
	    for (i in 1:N_fields){
	    	temp = eval(as.name(fields[i]))
			lty  = ltyvec[i]
			col  = colvec[i]
			points(fac*temp,1e-2*p_s,type="l",lwd=2,lty=lty,col=col)
	    } # field loop
	}  # mvals loop
	legend("topright",legend=fieldnames[1:N_fields],cex=cex, 
	      lty=ltyvec,col=colvec,lwd=2)

} #case loop
dev.off()
