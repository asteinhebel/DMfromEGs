#!/bin/bash
#Run likelihood stacking from SED outputs and interpret in DM space

#Assumes 	- directories for each target including data listed in file 'ft1.fits' 
#Inputs 	- Target name, black hole mass (1e8 solar masses) and distance (Mpc) 
#			_ Code pulls in computed SEDs from given targets
#Outputs 	-

# python schnittamn_likelihood_DM.py <target> <M_BH> <distance>

targets=("NGC4889" "NGC4649" "NGC1407")
m_bh=(208.0 47.3 46.5) #1e8 solar masses
distance=(102.0 16.5 29.0) #Mpc

#targets=("NGC4889")
#m_bh=(208.0) #1e8 solar masses
#distance=(102.0) #Mpc



for i in ${!targets[@]}; do
  python schnittman_likelihood_DM.py ${targets[$i]} ${m_bh[$i]} ${distance[$i]}
done

python plot_TSmaps_stack.py stackingTargets.txt