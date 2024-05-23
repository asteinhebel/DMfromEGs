#!/bin/bash
#Run likelihood stacking from SED outputs and interpret in DM space

#Assumes 	- directories for each target including data listed in file 'ft1.fits' 
#Inputs 	- Target name, black hole mass (1e8 solar masses) and distance (Mpc) 
#			_ Code pulls in computed SEDs from given targets
#Outputs 	- .npy file of DM likelihood values for each combination of DM mass/crossSection in <TARGET>/output/plots/dloglike
#			- likelihood plots in DM mass/crossSection space for every individual entry in <TARGET>/output/plots/
#			- likelihood plot of all entries in <TARGETS> stacked together in top level directory

targets=("dark4889")
m_bh=(1.0) #1e8 solar masses
distance=(99.03) #Mpc

for i in ${!targets[@]}; do
  python schnittman_likelihood_DM.py ${targets[$i]} ${m_bh[$i]} ${distance[$i]} -2 4 40 -30 -15 60
done

python plot_TSmaps_stack.py stackingTargets_dark.txt -2 4 40 -30 -15 60
