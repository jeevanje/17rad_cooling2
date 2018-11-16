Rtoolsdir = "~/Dropbox/Rtools/"
setwd("~/Dropbox/17rad_cooling2/")
source(paste(Rtoolsdir,"thermo_tools.R",sep=""))

D    = 1.5  # Diffusion parameter, as per RFM code for NQAD=1

#======#
# RFM  #
#======#

calc_cts2d = function(k,p_i,tabs_s,tau_i){
	nk        = dim(tau_i)[1]
	np        = dim(tau_i)[2]
	trans     = exp(-tau_i)
	dtrans_dp = (trans[ ,2:np]-trans[ ,1:(np-1)])/(rep(1,times=nk)%o%diff(p)) #s lev
	cts2d	  = -pi*outer(k,tabs_s,planck_k)*dtrans_dp  #s, W/m^2/Pa/m^-1, pppf units
	return(cts2d)
}

#===============#
# Simple models #
#===============#
load("./data/band_params.Rdata")  # includes ncoarse

calc_WVP0 = function(tabs,RH,Gamma,Ts=tabs[1]){
				# assume tabs on i levs, otherwise include Ts
				# Include D
				nz   = length(tabs)
		    	Tav  = 0.5*(tabs[1]+tabs[nz])
                WVP0 = D*RH*einf*Tav/L/Gamma
				return(WVP0)
}

calc_k1 = function(gas,p,tabs,Ts,Gamma,RH=0.75,q_co2=4.25e-4){
        if (gas == "h2o"){
            WVP0  = calc_WVP0(tabs,RH,Gamma,Ts)
            k1_rot  = k_rot + l_rot*(log(kappa_rot*WVP0)+log(ps/pref)+
                                 g/Rd/Gamma*log(tabs/Ts) - L/Rv/tabs) # s
            k1_vr   = k_vr  - l_vr*(log(kappa_vr*WVP0)+log(ps/pref)+
                                     g/Rd/Gamma*log(tabs/Ts) - L/Rv/tabs)
            k1_rot[k1_rot<k_rot] <- NA
            k1_vr[k1_vr>k_vr]    <- NA
			k1_list = list(k1_rot,k1_vr)
        } else if (gas=="co2"){
            k1_P = k_Q - l_Q*log(q_co2*p^2*kappa_Q/2/g/pref)
            k1_R = k_Q + l_Q*log(q_co2*p^2*kappa_Q/2/g/pref)
            k1_P[k1_P>k_Q] <- NA
            k1_R[k1_R<k_Q] <- NA
			k1_list = list(k1_P,k1_R)
        }
        return(k1_list)
}

calc_cts_verysimple = function(gas,p_s,tabs_s,Ts,Gamma,RH=0.75,q_co2=4.25e-4){
	np_s	         = length(p_s)
   	cts_simple_bands = array(dim=c(np_s,2))
	band_data_gas    = band_data[[gas]]
        if (gas == "h2o"){
            beta = 1 + L/Rv/tabs_s*Gamma*Rd/g
			#beta  = 1 + L/Rv/250*Gamma*Rd/g
			k1_list = calc_k1(gas,p_s,tabs_s,Ts,Gamma,RH,q_co2)
            k1_rot  = k1_list[[1]]
            k1_vr  = k1_list[[2]]
            cts_simple_bands[ ,1] = pi*planck_k(k1_rot,tabs_s)*l_rot*beta/p_s
            cts_simple_bands[ ,2] = pi*planck_k(k1_vr,tabs_s)*l_vr*beta/p_s
        } else if (gas=="co2"){
			beta    = 2
			k1_list = calc_k1(gas,p_s,tabs_s,Ts,Gamma,RH,q_co2)
            k1_P    = k1_list[[1]]
            k1_R    = k1_list[[2]]
            cts_simple_bands[ ,1] = pi*planck_k(k1_P,tabs_s)*l_Q*beta/p_s
            cts_simple_bands[ ,2] = pi*planck_k(k1_R,tabs_s)*l_Q*beta/p_s
        }
        cts_verysimple = apply(cts_simple_bands,1,sum,na.rm=TRUE)
        cts_verysimple[cts_verysimple==0] <- NA
		return(cts_verysimple)
}