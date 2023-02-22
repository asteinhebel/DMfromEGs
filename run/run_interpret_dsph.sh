#!/bin/bash
#Run likelihood stacking from SED outputs and interpret in DM space

#Assumes 	- directories for each target including data listed in file 'ft1.fits' 
#Inputs 	- Target name, J-factor (GeV^2/cm^5) and associated error. 
#			_ Code pulls in computed SEDs from given targets
#Outputs 	-

# python diMauro_likelihood_JfactorDM.py <target> <Jfactor> <JfactorError>


targets=("Draco" "Hercules" "SDG")
jfactor=(18.73 17.32 18.0)
error=(0.03 0.57 0.5)

for i in ${!targets[@]}; do
  python diMauro_likelihood_JfactorDM.py ${targets[$i]} ${jfactor[$i]} ${error[$i]}
done
