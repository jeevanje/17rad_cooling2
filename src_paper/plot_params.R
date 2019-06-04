source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")

# Labels
#gasnames  = c(expression(H[2]*O),expression(C[O]*2))
gasnames   = c("H2O","CO2")
Hk_units   = expression("("*K/d*a*y/c*m^{-1}*")")
Fklab      = expression(F[k]~~"("~W/m^2/c*m^{-1}~")") 
OLRklab    = expression(OLR[k]~~"("~W/m^2/c*m^{-1}~")") 
Hlab       = "H (K/day)"

# Factors
Nsec       = 3600*24   # sec/day
coo2d_fac  = -g/Cp*Nsec*1e2  # pppf to K/day/cm^-1 heating
coo1d_fac  = -g/Cp*Nsec      # pppf to K/day       heating
k_lim      = c(50,1450) # cm^-1  
p_lim      = c(1000,1) # hPa
coo2d_lim  = c(-0.014,0.014) # K/day/cm^-1

cex_leg    = 2.25

# plot functions
plot_coo2d = function(k,p_s,coo2d,main){
				# Assume s levs!
				cex  = 2.5
				np_s = length(p_s)
				pvec = np_s:1
				my.image.plot(1e-2*k,1e-2*p_s[pvec],coo2d_fac*coo2d[ ,pvec],
						xlim = k_lim,
						ylim = p_lim,
						zlim = coo2d_lim,
						xlab = klab,
						ylab = plab,
						main = main,
						cex.lab  = cex,
						cex.axis = cex,
						cex.main = cex,
						cex.legend = cex_leg,
						legend.width=1.75
						)
}

# Note: everything on s levs, since that is what's natural for Hk, and 
# want to show p1 curves on both Hk and tau plots and they should agree. 
plot_2dfield = function(k,z,z_fac,zlim,main,col=tim.colors(64),log=""){
	my.image.plot(1e-2*k,p_fac*p_s[p_vec],z_fac*z[ ,p_vec],
			xlim = k_lim,
			ylim = p_lim,
			zlim = zlim,
			#xlab = klab,
			ylab = p_lab,
			log  = log,
			main = main,
			col  = col,
			cex.lab  = cex,
			cex.axis = cex,
			cex.main = cex,
			legend.shrink = 0.9,
			legend.width  = 1.75,
			#legend.mar    = 5,
			cex.legend = cex_leg
		)
	title("")
	mtext(klab,side=1,line=4,outer=FALSE,cex=1.75)	

}

