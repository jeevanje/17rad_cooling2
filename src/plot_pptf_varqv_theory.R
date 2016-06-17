library(ncdf)
library(fields)

source("~/Dropbox/Rtools/thermo_tools.R")
source("~/Dropbox/Rtools/my_image_plot.R")

# Get data
Ts = 310
ncpath = "~/Dropbox/rad_cooling/data/prod_data_4_22_16/310k/data/verticalstats.nc"
nc = open.ncdf(ncpath)
z = get.var.ncdf(nc,"z")
nz = length(z)
time = get.var.ncdf(nc,"time")
nt = length(time)
nt_avg = 20
start = c(1,nt-nt_avg)
tabs = apply(get.var.ncdf(nc,start=start,"tabs"),1,mean)
p = apply(get.var.ncdf(nc,start=start,"p"),1,mean)
Ttp = 200
ktp = which.min(abs(tabs-Ttp))
ptp = p[ktp]
ps  = p[1]
pref_new = 5e4  # mid-troposphere
pref_old = 1e4
Tav = (Ttp+Ts)/2
gamma_av = g/Rd*log(Ttp/Ts)/log(ptp/ps)

# Simple theory
kappa0 = 50*pref_new/pref_old # m^2/kg, plus fudge factor for 75th percentile?
lk 	   =  6500 # m^-1
RHvec  = c(1, 0.8, 0.6, 0.4, 0.2)
N_rh   = length(RHvec) 
WVP0vec = RHvec*einf*Tav/L/gamma_av
k1_array = lk*outer(-L/Rv/tabs,log(kappa0*WVP0vec),"+")


#=========#
# Plots   #
#=========#
zvec=1:ktp
pptflim=c(0,5)
tabslim=c(Ts,200)
cex=1.5
colvec = rev(tim.colors(N_rh))
kvals = seq(100,999,length.out=100)     # cm^-1


pdf(file="../plots/pptf_varqv_theory.pdf", width = 11, height = 5)
par(mar=c(5,5,5,7),mfrow=c(1,2))

# Panel 1 -- ppptf profiles
plot(1,type="n",xlim=pptflim,ylim=tabslim,
        xlab=expression(-partialdiff[T]*F~~"("~W/m^2/K~")"),
        ylab="Temperature (K)",
        main = "LW Flux Divergence, theory",
        cex.main=cex,
        cex.lab=cex,
        cex.axis=cex)

for (i in 1:N_rh){
    pptf =  - pi*planck_k(tabs,k1_array[ ,i])*L/Rv/tabs^2*lk
    points(-pptf[zvec],tabs[zvec],type="l",lwd=2,col=colvec[i])
    }

legend("topright",legend = RHvec,title="RH",
        lty="solid",col=colvec,cex=1.25 )


# Panel 2 -- pptf_phase
pptf_function = function(T,k){
					pi*lk*planck_k(T,k)*L/Rv/T^2
					}
pptf_phase = t(outer(tabs,1e2*kvals,pptf_function))
my.image.plot(kvals,tabs[ktp:1], 
		  pptf_phase[ ,ktp:1],
		  zlim=range(pptf_phase[ ,ktp:1]),
      	  ylim = c(Ts,Ttp),
          main = expression(partialdiff[T]~F~~"("~W/m^2/K~")"),
          ylab = "Temperature (K)",
          xlab = expression("Wavenumber"~~"("~cm^-1~")"),
          cex.lab = cex,
          cex.axis = cex,
          cex.main = 1.5,
          cex.legend = cex)
for (i in 1:N_rh){
	points(1e-2*k1_array[zvec,i],tabs[zvec],type="l",lwd=2,col=colvec[i],lty="solid")
    }
#text(300,260,expression(k[1]*"("*T*")"),cex=1.25)


dev.off()
