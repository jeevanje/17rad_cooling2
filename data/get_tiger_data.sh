case=$1

# scripts
#scp nadirj@tiger.princeton.edu:~/rad_cooling2/R/plot_*.R ../src/

# data
#for case in co2_only_isotherm_atm co2_only_simple_atm h2o_no_cont_isotherm_atm h2o_only_no_cont
#do
   scp tigressdata:~/17rad_cooling2/rfm/${case}/${case}.nc ./
#done          

