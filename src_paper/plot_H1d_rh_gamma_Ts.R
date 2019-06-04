Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")

library(fields)
library(ncdf4)
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))
source(paste(Rtoolsdir,"rfm_tools.R",sep=""))
source(paste(Rtoolsdir,"plot_tools.R",sep=""))
source("./src_paper/plot_params.R")
source("./src_paper/cts_functions.R")
load("./data_paper/band_params.Rdata")  # includes ncoarse

# Data params
thermo_cases = c("h2o_noctm_Ts300_rh0.75_gamma7","h2o_noctm_Ts300_rh0.3_gamma7",
				 "h2o_noctm_Ts300_rh0.75_gamma5")
Tslist       = seq(270,300,by=10)
Ts_cases     = paste("h2o_noctm_Ts",Tslist,"_rh0.75_gamma7",sep="")
cases		 = c(thermo_cases,Ts_cases)
NTs 		 = length(Tslist)
Nthermo      = length(thermo_cases)

#=======#
# Data  #
#=======#

for (case in cases){
	ncpath  = paste("~/Dropbox/17rad_cooling2/data_paper/",case,".nc",sep="")
	nc      = nc_open(ncpath)
	k       = ncvar_get(nc,"k")
	dk      = k[2]-k[1]
	p       = ncvar_get(nc,"p")
	np      = length(p)
	z       = ncvar_get(nc,"z")
	tabs    = ncvar_get(nc,"tabs")
	Ts	    = tabs[1]
	Tstrat  = tabs[np]
	Tav     = 0.5*(Ts+Tstrat)
	q       = ncvar_get(nc,"q_h2o")
	Gamma   = -(diff(tabs)/diff(z))[1]
	RH      = (q/qsat(tabs,p))[1]	
	nk      = length(k)
	coo2d   = ncvar_get(nc,"coo")/abs(coo2d_fac)  # convert to SI 
	opt     = D*ncvar_get(nc,"opt")
	
	# coo, cts
	WVP      = calc_WVP0(tabs,RH,Gamma)*esat(tabs)  # i lev
	tau_simple = kappa_simple_h2o%o%(ps/pref*(tabs/Ts)^(g/Rd/Gamma)*WVP)
	p_s      = rfm_lni2s(p)
	z_s      = rfm_i2s(z)
	tabs_s   = rfm_i2s(tabs)
	coo2d_s  = rfm_i2s(coo2d)
	cts2d_s  = calc_cts2d(k,p,tabs_s,opt)
	cts_simple2d = calc_cts2d(k_fit_h2o,p,tabs_s,tau_simple)
	coo1d    = apply(coo2d_s,2,sum)*dk
	cts1d    = apply(cts2d_s,2,sum)*dk
	cts_simple1d = apply(cts_simple2d,2,sum)*dk_fit_h2o	
	cts_simple_bands = calc_cts_simple_bands("h2o",p_s,tabs_s,Ts,Gamma,RH)	
	cts_verysimple1d = apply(cts_simple_bands,1,sum,na.rm=TRUE)
	pvec     = 1:(min(which(cts_verysimple1d==0))-1)

	vars     = c("coo1d","cts1d","cts_verysimple1d","z_s",
				 "p_s","pvec","tabs_s","Ts","Gamma","RH")
	for (var in vars){
		assign(paste(var,case,sep="_"),eval(as.name(var)))	
	}
}


#==============#
# Plot prelims #
#==============#

# common params
coo1d_lim  = c(-3,0)   # K/day
p_fac	   = 1e-2
cex	       = 1.9
cex_leg    = 1.4
cex_leg2   = 1.5
lwd        = 2.5
pars      <- c('plt','usr')
Hlab       = "H (K/day)"
ltyvec     = c("solid","solid","solid")
colvec_Ts  = tim.colors(NTs)
colvec_thermo = c("black","blue","red")
legendvec_thermo= c("Base","RH=0.3",expression(Gamma==5~K/km))

plot_H_frame = function(main){
        plot(1,type="n", 
           xlim     = coo1d_lim, 
           ylim     = c(1000,0),
           xlab     = Hlab,
           ylab     = plab,
           main     = main,
           cex.lab  = cex,
           cex.axis = cex,
           cex.main = cex)      
}


#=======#
# Plot  #
#=======#
	
# PDF
file = "~/Dropbox/17rad_cooling2/plots/H1d_rh_gamma_Ts.pdf"
pdf(file,width=11,height=10,bg="white")
par(mfrow=c(2,2),mar=c(5,5,5,3))


# RFM var_thermo
plot_H_frame("(a)  RFM, variable RH, lapse")
for (n in 1:Nthermo){
	case = thermo_cases[n]
	p_s  = eval(as.name(paste("p_s",case,sep="_")))
	pvec = eval(as.name(paste("pvec",case,sep="_")))
	col  = colvec_thermo[n]
	coo1d= eval(as.name(paste("coo1d",case,sep="_")))
	cts_verysimple = eval(as.name(paste("cts_verysimple1d",case,sep="_")))
	
	coo1d[is.na(cts_verysimple)]<-NA	
    points(coo1d_fac*coo1d[pvec],p_fac*p_s[pvec],type="l",lwd=lwd,col=col)  
}
legend("bottomright",legend=legendvec_thermo,col=colvec_thermo,lwd=2,
		cex=cex_leg)

# Theory var_thermo 
plot_H_frame("(b)  SSM1D, variable RH, lapse")
for (n in 1:Nthermo){
	case = thermo_cases[n]
	p_s  = eval(as.name(paste("p_s",case,sep="_")))
	pvec = eval(as.name(paste("pvec",case,sep="_")))
	col  = colvec_thermo[n]
	cts_verysimple = eval(as.name(paste("cts_verysimple1d",case,sep="_")))
    points(coo1d_fac*cts_verysimple[pvec],p_fac*p_s[pvec],type="l",lwd=lwd,col=col)  
}
legend("bottomleft",legend=legendvec_thermo,col=colvec_thermo,lwd=2,
		cex=cex_leg)

# RFM varsst
plot_H_frame("(c)  RFM, variable Ts")
for (n in 1:NTs){
    Ts    = Tslist[n] 
    case  = Ts_cases[n]
    col   = colvec_Ts[n]
    tabs_s= eval(as.name(paste("tabs_s",case,sep="_")))
    p_s   = eval(as.name(paste("p_s",case,sep="_")))
	pvec = eval(as.name(paste("pvec",case,sep="_")))
    coo1d = eval(as.name(paste("coo1d",case,sep="_")))
    cts_verysimple1d = eval(as.name(paste("cts_verysimple1d",case,sep="_")))
    coo1d[is.na(cts_verysimple1d)]<-NA
    points(coo1d_fac*coo1d[pvec],p_fac*p_s[pvec],type="l",lwd=lwd,col=col)  
}
legend("topleft",legend=Tslist,lwd=2,cex=cex_leg2,col=colvec_Ts,
        title=expression(T[s]~~"(K)"))

# Theory varsst
plot_H_frame("(d)  SSM1D, variable Ts")
for (n in 1:NTs){
    Ts    = Tslist[n] 
    case  = Ts_cases[n]
    col   = colvec_Ts[n]
    tabs_s= eval(as.name(paste("tabs_s",case,sep="_")))
    p_s   = eval(as.name(paste("p_s",case,sep="_")))
	pvec = eval(as.name(paste("pvec",case,sep="_")))
    cts_verysimple1d = eval(as.name(paste("cts_verysimple1d",case,sep="_")))
    points(coo1d_fac*cts_verysimple1d[pvec],p_fac*p_s[pvec],type="l",lwd=lwd,col=col)  
}
legend("bottomleft",legend=Tslist,lwd=2,cex=cex_leg2,col=colvec_Ts,
        title=expression(T[s]~~"(K)"))


dev.off()	

