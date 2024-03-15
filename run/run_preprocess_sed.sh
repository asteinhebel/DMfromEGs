#!/bin/bash

#Run preprocessing script for EGs - creates parameter txt files that will be used to compute SEDs
#Compute and save SED with free background

#Assumes - SC/ltcube files saved in /Users/asteinhe/FermiLAT/BHinEGs_DM/data/
#		 -  data saved in /Users/asteinhe/FermiLAT/BHinEGs_DM/data/'+<TARGET>+'/ft1.fits' where each target has its own directory which contains the list of files titled ft1.fits 
#
#Outputs - Config file, run log, resulting fit parameters in directory named after the target
#		 - fits files from different fit models, fermipy files in <TARGET>/output
#		 - plots prefit and postfit of data/model and significance/excess in <TARGET>/output/plots

# python preprocess_likelihoodFitting.py <target> <ra> <dec>
# python get_sed.py <target>

#####EGs#####
targets=("NGC4889" "NGC4649" "NGC1407" "NGC3842" "NGC3091" "NGC1550")
ra=(195.034 190.917 55.0496 176.009 150.059 64.9081)
dec=(27.977 11.5526 -18.5804 19.9498 -19.6365 2.40998)

for i in ${!targets[@]}; do
  python preprocess_likelihoodFitting.py ${targets[$i]} ${ra[$i]} ${dec[$i]}
  python get_sed.py ${targets[$i]}
done