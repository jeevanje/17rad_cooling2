rm(list=ls())  # clean up workspace before saving
load("~/Dropbox/17rad_cooling2/data_paper/band_data.Rdata")  # includes ncoarse

# h20
band_data_gas = band_data[['h2o']]
rot_data      = band_data_gas[[1]]
vr_data       = band_data_gas[[2]]
k_fit_rot     = rot_data$k_fit
k_rot	      = min(k_fit_rot)
kappa_rot     = max(rot_data$kappa_fit)
l_rot	      = rot_data$lk
k_fit_vr      = vr_data$k_fit
k_vr	      = max(k_fit_vr)
kappa_vr      = max(vr_data$kappa_fit)
l_vr	      = vr_data$lk
k_fit_h2o     = c(k_fit_rot,k_fit_vr[-1])
nk_fit_h2o    = length(k_fit_h2o)
dk_fit_h2o    = diff(k_fit_h2o)[1]
kappa_simple_h2o  = c(kappa_rot*exp(-abs(k_fit_rot -k_rot)/l_rot),
                   kappa_vr*exp(-abs(k_fit_vr[-1] -k_vr)/l_vr))

# co2
band_data_gas = band_data[['co2']]
P_data        = band_data_gas[[1]]
R_data        = band_data_gas[[2]]
k_fit_P       = P_data$k_fit
k_P  	      = max(k_fit_P)
kappa_P       = max(P_data$kappa_fit)
l_P  	      = P_data$lk
k_fit_R       = R_data$k_fit
k_R	          = min(k_fit_R)
kappa_R       = max(R_data$kappa_fit)
l_R  	      = R_data$lk
k_fit_co2     = c(k_fit_P,k_fit_R[-1])
nk_fit_co2    = length(k_fit_co2)
dk_fit_co2    = diff(k_fit_co2)[1]

k_Q  	      = mean(c(k_P,k_R))
l_Q	          = mean(c(l_P,l_R))		
kappa_Q	      = exp(mean(c(log(kappa_P),log(kappa_R))))
kappa_simple_co2  = kappa_Q*exp(-abs(k_fit_co2-k_Q)/l_Q)

save(list=ls(),file="~/Dropbox/17rad_cooling2/data_paper/band_params.Rdata")

