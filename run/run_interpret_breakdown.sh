#!/bin/bash
#Run likelihood stacking from SED outputs and interpret in DM space

#Assumes 	- directories for each target including data listed in file 'ft1.fits' 
#Inputs 	- Target name, black hole mass (1e8 solar masses) and distance (Mpc) 
#			_ Code pulls in computed SEDs from given targets
#Outputs 	- .npy file of DM likelihood values for each combination of DM mass/crossSection in <TARGET>/output/plots/dloglike
#			- likelihood plots in DM mass/crossSection space for every individual entry in <TARGET>/output/plots/
#			- likelihood plot of all entries in <TARGETS> stacked together in top level directory


targets=("NGC4889" "NGC4649" "NGC1407" "NGC3842" "NGC3091" "NGC1550" "NGC1600")
m_bh=(208.0 47.3 46.5 90.0 37.2 38.7 170.0) #1e8 solar masses
#distanceKH=(102.0 16.5 29.0 92.2 53.02 52.5) #Mpc
distance=(99.03 21.14 24.44 96.78 63.55 53.58 68.15) #Mpc


for i in ${!targets[@]}; do
  python schnittman_likelihood_DM.py ${targets[$i]} ${m_bh[$i]} ${distance[$i]} -2 4 40 -30 -15 60
done

python plot_TSmaps_stack.py stackingTargets.txt -2 4 40 -30 -15 60