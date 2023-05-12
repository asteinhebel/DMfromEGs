#!/bin/bash
#Run likelihood stacking from SED outputs and interpret in DM space

#Assumes 	- directories for each target including data listed in file 'ft1.fits' 
#Inputs 	- Target name, J-factor (GeV^2/cm^5) and associated error. 
#			_ Code pulls in computed SEDs from given targets
#Outputs 	-

# python diMauro_likelihood_JfactorDM.py <target> <Jfactor> <JfactorError>

targets=("NGC4889" "NGC4649" "NGC1407")
jfactor=(20 20 20) #exponent of it
#jfactor=(17 17 17) #exponent of it
error=(1 1 1) #log10J_error


for i in ${!targets[@]}; do
  python diMauro_likelihood_JfactorDM.py ${targets[$i]} ${jfactor[$i]} ${error[$i]}
done

