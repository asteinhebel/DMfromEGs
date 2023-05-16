#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from scipy import interpolate
from iminuit import Minuit
from scipy.optimize import fmin
import os,sys
import importlib
import copy
import yaml, shutil

############################################################################################################
# HELPER FUNCTIONS
############################################################################################################
def flux_DM_model(DMprop):
    
    #Compute dNdE DM from model
    DMmass = DMprop[0]
    sigmav = DMprop[1]  
    mBH = DMprop[2]
    dist = DMprop[3]
    
    sig26 = sigmav * 1e26 #scale to get sigma ^-26
    
    return (2.e-20)*(mBH**(5./2.))*sig26*pow(1./DMmass,2.)**pow(1./dist,2.)

def funcLogLike(DMprop,table_norm,table_loglike, eref_vec):
    #print(f"RUN AT [DMmass sigmav mBH dist] = {DMprop}")
    
    LogLike = 0.
    fluxval = flux_DM_model(DMprop)       
    #print(f"FLUX VAL = {fluxval}")
    
    #calculate individually in each energy bin - assume full flux is in every bin (so over-estimating, especially in high-energy bins)
    for t in range(len(table_norm)):
        #approximate function SED loglike = y = f(x) = f(table_norm)
        loglike_interp = interpolate.interp1d(table_norm[t],table_loglike[t], bounds_error=False, fill_value='extrapolate')
     
        #print(f"dloglike interp(fluxval) = {loglike_interp(fluxval)}")   
        
        """   
        #visualize interpolation in each bin
        xrn =  np.sort(table_norm[t])
        plt.plot(table_norm[t], table_loglike[t], 'o', xrn, loglike_interp(xrn), '-', fluxval, loglike_interp(fluxval), 'o')
        plt.xlabel(f'measured flux, energy bin {t}')
        plt.ylabel(f'dloglikelihood, energy bin {t}')
        plt.show()
        """
        
        LogLike += loglike_interp(fluxval)
        #print(f"{eref_vec[t]} GeV bin, {fluxval} /cm^2/s flux, {-loglike_interp(fluxval)} dloglike")
        
    return LogLike

def makeDir(dirpath):
    if(os.path.isdir(dirpath)==True):
        print(f"Path to {dirpath} exists")
        shutil.rmtree(dirpath)
        os.system('mkdir -p %s' %dirpath) 
    else:
        print(f"path to {dirpath} does not exist, creating...")
        os.system('mkdir -p %s' %dirpath)

############################################################################################################
# MAIN
############################################################################################################
def main(cmd_line):

	'''Here we can set variable for path definitions
		In addition to all desired parameters,
		double check you have the correct LTcube, extdir
	'''

	print("*****running!!!*****")
    
	srcname = cmd_line[1]
	##Jfactor = 10**float(cmd_line[2]) #J-factor in units of GeV^2/cm^5
	##sigmaJ = float(cmd_line[3]) #log10J_error
	m_BH = float(cmd_line[2]) #black hole mass, units of 10e8 solar masses
	distance = float(cmd_line[3]) #distance to galacy, units of Mpc


	sed_suffix = 'free_BG_sed'
	homedir = "/Users/asteinhe/FermiLAT/BHinEGs_DM/run/"

	sedfits =  f"{homedir}{srcname}/output/{srcname}_{sed_suffix}.fits"
	save_array_path = f"{homedir}{srcname}/output/dloglike/test/"
	makeDir(save_array_path)

	#print("Running for {}, J: {:.2e}, sigmaJ: {}".format(srcname,Jfactor, sigmaJ))
	print("Running for {}".format(srcname))

	sed = fits.open(sedfits)[1].data
	table_refdnde = copy.deepcopy(sed['ref_dnde']) #Differential flux of the reference model evaluated at the bin center (cm−2 s−1 MeV−1)
	table_normscan = copy.deepcopy(sed['norm_scan']) #Array of NxM normalization values for the profile likelihood scan in N energy bins and M scan points. A row-wise multiplication with any of ref columns can be used to convert this matrix to the respective unit.
	#convert table_normscan to units of cm^-2 s^-1 MeV^-1
	for t in range(len(table_refdnde)):
		table_normscan[t] = sed['norm_scan'][t,:]*table_refdnde[t]

	table_loglikescan = copy.deepcopy(sed['dloglike_scan']) #Array of NxM delta-loglikelihood values for the profile likelihood scan in N energy bins and M scan points.
	
	# Convert MeV to GeV
	print("Define eref_vec in GeV...")
	eref_vec = copy.deepcopy(sed['e_ref'])/1.e3 #upper edges of SED energy bins in GeV

	#-------------------------------------------
	# Producing Loglike profiling with sigmav.
	#-------------------------------------------

	# create results file
	print("create_vectors")
	mass_vec =  np.logspace(0,4,40) #GeV
	sigmav_vec = np.logspace(-28,-22,60)#cm^3/s
	LogLike_vec = np.zeros(shape=(len(mass_vec),len(sigmav_vec)))

	for u in range(len(mass_vec)):
		for t in range(len(sigmav_vec)):
			DMprop = [mass_vec[u],sigmav_vec[t], m_BH, distance]
			LogLike_vec[u,t] = funcLogLike(DMprop,table_normscan,table_loglikescan, eref_vec=eref_vec)
			#print('{:.2e} {:.2f}, {}'.format(sigmav_vec[t], mass_vec[u], np.round(LogLike_vec[u],2)))

	np.save(save_array_path+'{}_freeBG_dlike.npy'.format(srcname), -LogLike_vec)

################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":
    main(sys.argv)

