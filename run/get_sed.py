import os,sys 
import matplotlib.pyplot as plt
import numpy as np
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter
import yaml


def main(cmd_line):

    '''
    	Calculates the SED for the src given by the command line argument
    '''
    


    print("*****running!!!*****")
    srcname = cmd_line[1]
    preprocessing_path = "/Users/asteinhe/FermiLAT/BHinEGs_DM/run/%s/"%(srcname)

    param_file = os.popen("ls "+preprocessing_path+"*_Param.txt").read()
    param_name = (param_file.split('/')[-1]).split('_Param.txt')[0]
    
    if srcname!=param_name:
    	print( srcname + ' ' + param_name)
    else:
    	print(srcname)

    
    print( "change dir...")
    os.chdir(preprocessing_path)

    print( "initialize GTA...")
    gta = GTAnalysis(preprocessing_path+srcname + ".yaml", logging={"verbosity":3})

    print( "change dir...")
    os.chdir("output")
    gta.load_roi("fit_model_3.npy")

    print( "Starting fit...")
    gta.optimize()

    print( "starting sed...")

    
    #gta.sed(param_name,write_npy=True, write_fits=True, free_background=True,  outfile=srcname+"_free_BG_sed.fits")
    gta.sed(param_name, free_background=True,  outfile=srcname+"_free_BG_sed.fits", make_plots=True)

    print( "Job complete!")
    

################################################
if __name__=="__main__":
    main(sys.argv)




