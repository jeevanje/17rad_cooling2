load("~/Dropbox/17rad_cooling2/data/band_data.Rdata")  # includes ncoarse

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

k_Q  	      = mean(c(k_P,k_R))
l_Q	      = mean(c(l_P,l_R))		
kappa_Q	      = exp(mean(c(log(kappa_P),log(kappa_R))))

save(list=ls(),file="~/Dropbox/17rad_cooling2/data/band_params.Rdata")

