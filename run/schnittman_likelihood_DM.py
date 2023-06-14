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
def diffEn(m,e):
	frac = (e/m)**2
	flux = frac * (e/m) * np.exp(-frac)
	return flux    

def full_diffPhotSpectrum(m, mBH, dist): 
    #Calculate dE/dN of DM spectrum, units of ph/cm^2/s/GeV (after multiplication with DM sigma), return as array of points that span LAT sensitivity 
    norm = (2.e-20)*pow(mBH,(5./2.))*pow(1./m,2.)*pow(1./dist,2.) * 2 / m
    spectrum = np.array([norm*diffEn(m, en) for en in obsEnergies])
    return spectrum

def flux_DM_model(DMprop, spec, eref, emin):
    #Compute dNdE DM from model
    DMmass = DMprop[0]
    sigmav = DMprop[1]  
    mBH = DMprop[2]
    dist = DMprop[3]
    
    #Finish normalizing photon spectrum
    sig26 = sigmav * 1e26 #scale to get sigma ^-26
    
    fullSpec = np.array([s*sig26 for s in spec]) #= dN/dE [ph/cm^2/s/GeV]
    interp_spectrum = interpolate.InterpolatedUnivariateSpline(obsEnergies, fullSpec, k=1, ext=0) #linear interpolation because k=1, ext=0=>extrapolate. Identical to interp1d
    
    #Integrate photon flux in each SED bin from min-max value (emin to eref)
    fluxInBins=[]
    for sedbin in range(len(eref)):
        #define finer mesh of points within SED bin to compute integral
        xx = np.linspace(emin[sedbin], eref[sedbin], 200)
        #integrate interpolated flux over this range
        fluxInBins.append(interp_spectrum.integral(xx[0], xx[-1]))

    fluxInBins = np.array(fluxInBins)
    
    """
    #Confirm spectrum shape/scale is as expected
    print(f"DM MASS - {DMmass} GeV")
    print(f"DM SIGMA - {sig26} e26 cm^3/s")
    #print(f"Flux at SED edges - {fluxInBins}")
    if DMmass>90 and DMmass<120 and sig26>0.8 and sig26<2:
        plt.plot(eref, fluxInBins, '-o')
        plt.xlabel("GeV")
        plt.ylabel("ph/cm^2/s")
        plt.xscale('log')
        plt.yscale('log')
        plt.show()
    """
    
    return fluxInBins

def funcLogLike(DMprop, spec, table_norm, table_loglike, eref_vec, emin_vec):
    #print(f"RUN AT [DMmass sigmav mBH dist] = {DMprop}")
    
    LogLike = 0.
    fluxval = flux_DM_model(DMprop, spec, eref_vec, emin_vec) 
    #array the same length as number of SED bins, each value has units of ph/cm^2/s     
    #print(f"FLUX VAL = {fluxval}") 

    #calculate dloglike of DM flux value individually in each energy bin from interpolated data
    for t in range(len(table_norm)):
        #approximate function SED loglike = y = f(x) = f(table_norm)
        loglike_interp = interpolate.interp1d(table_norm[t],table_loglike[t], bounds_error=False, fill_value='extrapolate')

        print(f"dloglike interp(fluxval) = {loglike_interp(fluxval)}")   
        
        """
        #FOR DEBUGGING 
        #visualize interpolation in each bin
        print(f"bin {t}, DM flux = {fluxval[t]} with dloglike={loglike_interp(fluxval[t])}")
        xrn =  np.sort(table_norm[t])
        xrn_longer = np.logspace(-25,-12, 100)
        plt.plot(table_norm[t], table_loglike[t], 'o', xrn_longer, loglike_interp(xrn_longer), '-', fluxval[t], loglike_interp(fluxval[t]), 'o')
        plt.xlabel(f'measured flux [ph/cm^2/s/MeV], energy bin {t}')
        plt.ylabel(f'dloglikelihood, energy bin {t}')
        plt.xscale('log')
        plt.ylim((-25, 2))
        plt.show()
        plt.clf()
        """

        LogLike += loglike_interp(fluxval[t])
        
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
	#convert table_normscan to units of cm^-2 s^-1 GeV^-1
	for t in range(len(table_refdnde)):
		table_normscan[t] = sed['norm_scan'][t,:]*table_refdnde[t]/1.e3

	table_loglikescan = copy.deepcopy(sed['dloglike_scan']) #Array of NxM delta-loglikelihood values for the profile likelihood scan in N energy bins and M scan points.
	
	# Convert MeV to GeV
	print("Define eref_vec in GeV...")
	eref_vec = copy.deepcopy(sed['e_ref'])/1.e3 #upper edges of SED energy bins in GeV
	emin_vec = copy.deepcopy(sed['e_min'])/1.e3 #lower edges of SED energy bins in GeV

	#-------------------------------------------
	# Producing Loglike profiling with sigmav.
	#-------------------------------------------

	# create results file
	print("create_vectors")
	#mass_vec =  np.logspace(0,4,40) #GeV
	mass_vec =  np.logspace(-1,4,40) #GeV
	sigmav_vec = np.logspace(-28,-22,60)#cm^3/s
	LogLike_vec = np.zeros(shape=(len(mass_vec),len(sigmav_vec)))

	for u in range(len(mass_vec)):
	    #plot energy spectrum with all normalization except XS
	    photSpectrum = full_diffPhotSpectrum(mass_vec[u], m_BH, distance)
	    for t in range(len(sigmav_vec)):
	        DMprop = [mass_vec[u],sigmav_vec[t], m_BH, distance]
	        LogLike_vec[u,t] = funcLogLike(DMprop,photSpectrum,table_normscan,table_loglikescan,eref_vec,emin_vec)
			#print('{:.2e} {:.2f}, {}'.format(sigmav_vec[t], mass_vec[u], np.round(LogLike_vec[u],2)))

	np.save(save_array_path+'{}_freeBG_dlike.npy'.format(srcname), LogLike_vec)

################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":

    #LAT data range
    obsEnergies = np.logspace(0,3,200) #GeV
    #obsEnergies = [1, 10]
    main(sys.argv)

