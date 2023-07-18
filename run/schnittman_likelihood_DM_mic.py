
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
print("create_vectors")
obsEnergies = np.logspace(-1,3,200) #GeV
mass_vec =  np.logspace(-1,3,10) #GeV
sigmav_vec = np.logspace(-28,-22,10) #cm^3/s


def diffEn(m,e):
    frac = (e/m)**2
    flux = frac * (e/m) * np.exp(-frac)
    return flux

def full_diffPhotSpectrum(m, sv, mBH, dist):
    #Calculate dE/dN of DM spectrum, units of ph/cm^2/s/GeV (after multiplication with DM sigma), return as array of points that span LAT sensitivity
    norm = 2.e-20 * pow(mBH, (5./2.)) * (sv/1e-26) * pow(1./m, 2.) * pow(1./dist,2.) * (2 / m) # ph/cm2/s/GeV
    spectrum = np.array([norm * diffEn(m, en) for en in obsEnergies])
    
    return spectrum
    
def generate_flux_DM_model_file(DMprop, spec, folder):

    dict = {'Energy[MeV]'           : obsEnergies*1000,
            'DiffFlux[ph/cm2/s/GeV]': spec
            }
            
    df = pd.DataFrame.from_dict(dict,)
    df.to_csv('%s/Spec_DM_m%.1f_sv%i.txt'%(folder, DMprop[0], -np.log10(DMprop[1])),  header=False, index=False)
    return 0


def create_DMspectra_files(cmd_line):
    '''
    '''
    print("Creating the files")
    srcname = cmd_line[1]
    m_BH = float(cmd_line[2])
    distance = float(cmd_line[3])
    
    sed_suffix = 'free_BG_sed'
    homedir = "/Users/mnegro/MyDocuments/_FERMI/DM_EllipticalGalaxies/DMfromEGs/run/"
    
    src_folder = f"{homedir}{srcname}/"
    sedfits =  f"{homedir}{srcname}/output/{srcname}_{sed_suffix}.fits"
    
    for u in range(len(mass_vec)):
    
        for t in range(len(sigmav_vec)):
            print('DM mass, DM xsection: ', mass_vec[u], sigmav_vec[t])
            photSpectrum = full_diffPhotSpectrum(mass_vec[u], sigmav_vec[t], m_BH, distance)
            DMprop = [mass_vec[u], sigmav_vec[t], m_BH, distance]
            generate_flux_DM_model_file(DMprop, photSpectrum, src_folder)


################################################################################
################################################################################
################################################################################
if __name__=="__main__":

    #------- Mik test
    targets = ["NGC4889", "NGC4649", "NGC1407"]
    ra      = [195.034, 190.917, 55.0496    ]
    dec     = [27.977, 11.5526, -18.5804    ]
    m_bh    = [208.0, 47.3, 46.5] #1e8 solar masses
    distance= [102.0, 16.5, 29.0] #Mpc
    
    create_DMspectra_files(['placeholder', targets[0], m_bh[0], distance[0]])


