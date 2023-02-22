A typical run of these scripts for a DM analysis usually involves 
i) running the preprocessing step (preprocess_JLA.py) 
ii) computing the sed (get_sed.py) 
iii) running the diMauro_likelihood.py script (diMauro_likliehood.py) to convert the the delta-log-likelihood from flux-energy space (i.e. the SED) into mass-sigmav space. 

Some details on the included files:
project_config.yaml - defines the relevant directory paths, data sources and LT Cube path.
preprocess_JLA.py - the script for our preprocessing step. This does the Fermi joint likelihood analysis and does not have anything to do with DM. Required inputs are some source name, RA (deg), and Dec (deg)
get_sed.py - this loads the results of the preprocessing step and calculates the SED
diMauro_likliehood.py - This is primarily Mattia’s and uses the results of the SED calculation to produce a delta-log-likelihood profile in mass-sigmav space. This file takes as input the target name as well as the J-factor and J-factor uncertainty (which it sounds that you ultimately wont need). It also requires the DM particle spectrum, which for our analysis we use the ‘AtProductionNoEW_gammas.dat’ file that Ive also included.

Everything is written in python 3 and uses fermipy v1.2. All the path directories are setup for our computing cluster, so those you’ll need to make changes to work with your system. The exctractcirellitable() function in the diMauro_likliehood.py script provides the gamma-ray spectrum from WIMP DM annihilation, so youll also need to change this to work with your DM models.
