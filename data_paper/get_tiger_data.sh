case=$1

# scripts
#scp nadirj@tiger.princeton.edu:~/rad_cooling2/R/plot_*.R ../src/

# data
#for case in h2o_noctm_Ts300_rh0.75_gamma5 h2o_noctm_Ts290_rh0.75_gamma7 h2o_noctm_Ts280_rh0.75_gamma7 h2o_noctm_Ts270_rh0.75_gamma7 ; 
#do
   scp tigressdata:~/17rad_cooling2/rfm/${case}/${case}.nc ./
#done          

