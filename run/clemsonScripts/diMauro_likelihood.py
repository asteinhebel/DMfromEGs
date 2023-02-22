#!/usr/bin/env python
# coding: utf-8

# In[50]:


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
import yaml


# Function Defs
def setplot():
    '''Defines standard plot parameters'''
    
    fig = plt.figure(figsize=(3.5, 3), dpi=150)
    ax = plt.gca()
    ax.tick_params(direction='in', which='both', top=True, right=True)
    return ax, fig




def exctractcirellitable(DMmass,DMchannel,particle,EWcorr):
    '''
    This function returns the energy spectrum in 1/GeV for the particle production from DM annihilation.
    DMmass: dark matter mass in GeV
    DMchannel: dark matter annihilation channel
        #e 4 (2), mu 7 (3), tau 10 (4), bb 13 (7), tt 14 (8), WW 17 (9), ZZ 20 (10), gamma 22 (12), h 23 (13)
    particle: particle produced from the DM annihilation ('gammas' or 'positrons')
    EWcorr: electroweak corrections ('Yes' or 'No')
    '''
    
    #with EW corrections
    #mDM      Log[10,x]   eL         eR         e          \[Mu]L     \[Mu]R     \[Mu]      \[Tau]L    \[Tau]R    \[Tau]     q            c            b            t            WL          WT          W           ZL          ZT          Z           g            \[Gamma]    h           \[Nu]e     \[Nu]\[Mu]   \[Nu]\[Tau]   V->e       V->\[Mu]   V->\[Tau]
    
    #NO EW Corr.
    #   mDM      Log[10,x]   e           \[Mu]      \[Tau]     q           c           b           t           W          Z          g           \[Gamma]   h
    
    if EWcorr=='Yes':
        listenergies = 179
        energy_vec = np.arange(-8.9,0.05,0.05)
    elif EWcorr=='No':
        listenergies = 180
        energy_vec = np.arange(-8.95,0.05,0.05)
    else:
        print('Error Wrong value for EWcorr, Yes or No')
    
    energy = np.zeros(listenergies)
    fluxDM = np.zeros(listenergies)
    
    table = np.loadtxt('/home/armcdan/dSphs/Alt_stacking/particle_data/AtProduction%sEW_%s.dat'%(EWcorr,particle), skiprows=1)
    
    massvec = []
    for t in range(len(table)):
        if t%listenergies == 0:
            massvec.append(table[t,0])
    massvec = np.array(massvec)
    
    flux = []
    for t in range(len(table)):
        flux.append(table[t,DMchannel])
    
    f = interpolate.interp2d(massvec, energy_vec, flux, kind='linear')

    for t in range(len(energy_vec)):
        fluxDM[t] = f(DMmass,energy_vec[t])
    
    return np.power(10.,energy_vec)*DMmass,fluxDM/(np.log(10.)*np.power(10.,energy_vec)*DMmass)

def func_interpolate(varval,variablevec,funcvec):
    
    result = 0.
    if varval<variablevec[0]:
        result = 0.
    elif varval>variablevec[len(variablevec)-1]:
        result = 0.
    else:
        Log10E_bin = np.log10(variablevec[1])-np.log10(variablevec[0])
        nbin = ( np.log10(varval)-np.log10(variablevec[0]) )/Log10E_bin
        binval = int(nbin)
        #print(fluxDM[binval],fluxDM[binval+1],EnergyDM[binval],EnergyDM[binval+1])
        result = pow(10.,  np.log10(funcvec[binval]) + (np.log10(funcvec[binval+1])-np.log10(funcvec[binval]))*(np.log10(varval)-np.log10(variablevec[binval]))/(np.log10(variablevec[binval+1])-np.log10(variablevec[binval]))  )
        #print(Energy,EnergyDM[binval],EnergyDM[binval+1],fluxDM[binval],fluxDM[binval+1],result)
    #print('Interpolate',varval,variablevec[binval],variablevec[binval+1],funcvec[binval],funcvec[binval+1],result)
    return result

def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
    '''
    Adapted and modifed to get the unknowns for defining a parabola:
    https://urldefense.com/v3/__http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points__;!!PTd7Sdtyuw!UQyjBqduMfb3gwgoTAUvWDfRQarQOwKMc8x3FXyhvPrN7xITpDrqjHfbo4qSiQuOyoiN9DeUH81_c-kTwwYM-Tjvuw$  
    '''

    denom = (x1-x2) * (x1-x3) * (x2-x3);
    A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom;
    B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom;
    C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom;
    return A,B,C


def interpolate_func(sed_scan,dloglike_scan,ris):
    if ris<sed_scan[0]:
        return dloglike_scan[0]+(dloglike_scan[1]-dloglike_scan[0])/(sed_scan[1]-sed_scan[0])*(ris-sed_scan[0])
    elif ris>sed_scan[len(sed_scan)-1]:
        return dloglike_scan[len(sed_scan)-2]+(dloglike_scan[len(sed_scan)-1]-dloglike_scan[len(sed_scan)-2])/(sed_scan[len(sed_scan)-1]-sed_scan[len(sed_scan)-2])*(ris-sed_scan[len(sed_scan)-2])
    else:
        f = interp1d(sed_scan, dloglike_scan)
        #denergy = np.log10()
        #print sed_scan,dloglike_scan,ris
        #print f(ris)
        return f(ris)

#sed_scan[t],loglike_table[t][:],sed_stacking_best[t]
def findLogLprofile_func(loglike_table,sed_scan):
    a,b,c=calc_parabola_vertex(sed_scan[0], loglike_table[0], sed_scan[1], loglike_table[1], sed_scan[2], loglike_table[2])
    sed_best = -b/(2.*a)
    LogL_best = a*sed_best*sed_best+b*sed_best+c
    deltaLogL = LogL_best-1.
    deltaparabola = b*b - 4.*a*(c-deltaLogL)
    error_down = (-b-sqrt(deltaparabola))/(2.*a)
    error_up = (-b+sqrt(deltaparabola))/(2.*a)
    #print sed_best,LogL_best,error_up,error_down,a,b,c,sed_scan[0], loglike_table[0], sed_scan[1], loglike_table[1], sed_scan[2], loglike_table[2]
    return sed_best,LogL_best,error_up,error_down,a,b,c

def flux_DM_prompt(Energy,DMprop,Jfact):
    
    #Precompute dNdE DM
    DMmass = DMprop[0]
    sigmav = DMprop[1]
    DMchannel = DMprop[2]
    EWcorr = DMprop[3]    
    EnergyDM,fluxDM = exctractcirellitable(DMmass,DMchannel,'gammas',EWcorr)
    for t in range(len(EnergyDM)):
        if fluxDM[t]>0.:
            fluxDM[t]=fluxDM[t]
        else:
            fluxDM[t]=1e-30
    
    dNdE = func_interpolate(Energy,EnergyDM,fluxDM)

    return 0.5*(1./(4.*np.pi))*pow(1./DMmass,2.)*Jfact*sigmav*dNdE

#####

def interpolate_Loglike(fluxval,table_norm,table_loglike):
    
    #print(table_norm,table_loglike)
    
    if fluxval>table_norm[len(table_norm)-1]:
        return table_loglikescan[len(table_norm)-1]
    else:
        f = interpolate.interp1d(table_norm,table_loglike, kind='linear')
        return f(fluxval)
    
def funcLogLike_old(DMprop,table_norm,table_loglike, eref_vec):
    
    LogLike = 0.
    for t in range(len(table_norm)):
        #tstart_0 = timer()
        fluxval = flux_DM_prompt(eref_vec[t],DMprop,Jfactor)/1.e3
        if fluxval<table_norm[t][len(table_norm[t])-1]:
            LogLike = LogLike + interpolate_Loglike(fluxval,table_norm[t],table_loglike[t])
        else:
            LogLike = LogLike + table_loglike[t][len(table_loglike[t])-1]
        #print(fluxval,table_norm[t],table_loglike[t],LogLike)
        #print('Time',abs(tstart_0-timer()))
    return -LogLike


def funcLogLike(DMprop,table_norm,table_loglike, eref_vec, Jfactor):
    
    LogLike = 0.
    for t in range(len(table_norm)):

        fluxval = flux_DM_prompt(eref_vec[t],DMprop,Jfactor)/1.e3
        loglike_interp = interpolate.interp1d(table_norm[t],table_loglike[t], bounds_error=False, fill_value='extrapolate')
        
        LogLike += loglike_interp(fluxval)
        
    return -LogLike

def funcLogLike_profJ_old(DMprop,table_norm,table_loglike, eref_vec):
    
    fluxvec = np.zeros(len(table_normscan))
    for t in range(len(table_normscan)):
        fluxvec[t] = flux_DM_prompt(eref_vec[t],DMprop,Jfactor)/1.e3
        
    def f(par0):
        
        Loglike = 0
        LogLikeJ = np.log10(np.exp(-pow((np.log10(Jfactor*par0)-np.log10(Jfactor)),2.)/(2.*sigmaJ*sigmaJ)))
        
        for t in range(len(table_normscan)):
            fluxval = fluxvec[t]*par0
            
            if fluxval<table_norm[t][len(table_norm[t])-1]:
                Loglike = Loglike + interpolate_Loglike(fluxval,table_norm[t],table_loglike[t])
            else:
                Loglike = Loglike + table_loglike[t][len(table_loglike[t])-1]
            #Loglike = Loglike + interpolate_Loglike(fluxval,table_norm[t],table_loglike[t])
        
        #print(par0,-Loglike,-LogLikeJ,-Loglike-LogLikeJ)
        return -Loglike-LogLikeJ

    #m=Minuit(f, par0=1., error_par0=1e-4, limit_par0=(0.10,1.00), print_level=1,errordef=1)
    m=Minuit(f, par0=1.)
    m.limits = (0.01,10)
    m.migrad()
    
    Jren_bestfit = m.values["par0"]
    Jren_error = m.errors["par0"]
    
    return m.fval,Jren_bestfit,Jren_error


def funcLogLike_profJ(DMprop,table_norm,table_loglike, eref_vec, Jfactor, sigmaJ):
    
    fluxvec = np.zeros(len(table_norm))
    for t in range(len(table_norm)):
        fluxvec[t] = flux_DM_prompt(eref_vec[t],DMprop,Jfactor)/1.e3
        
    def f(par0):
        
        Loglike = 0
        LogLikeJ = np.log10(np.exp(-pow((np.log10(Jfactor*par0)-np.log10(Jfactor)),2.)/(2.*sigmaJ*sigmaJ)))
        
        for t in range(len(table_norm)):
            fluxval = fluxvec[t]*par0
            
            if fluxval<table_norm[t][len(table_norm[t])-1]:
                Loglike = Loglike + interpolate_Loglike(fluxval,table_norm[t],table_loglike[t])
            else:
                Loglike = Loglike + table_loglike[t][len(table_loglike[t])-1]
            #Loglike = Loglike + interpolate_Loglike(fluxval,table_norm[t],table_loglike[t])
        
        #print(par0,-Loglike,-LogLikeJ,-Loglike-LogLikeJ)
        return -Loglike-LogLikeJ

    #m=Minuit(f, par0=1., error_par0=1e-4, limit_par0=(0.10,1.00), print_level=1,errordef=1)
    m=Minuit(f, par0=1.)
    m.limits = (0.1,10) #changed from 0.01,10
    m.errors = 1e-4 #added
    m.tol=1e-4 #added
    m.errordef=1 #added
    m.migrad()
    
    Jren_bestfit = m.values["par0"]
    Jren_error = m.errors["par0"]
    
    return m.fval,Jren_bestfit,Jren_error




def main(cmd_line):

    '''Here we can set variable for path definitions
       In addition to all desired parameters,
       double check you have the correct LTcube, extdir
    '''
    with open("project_config.yaml", "r") as stream:
        config_dict = yaml.safe_load(stream)

    cuser = config_dict['cuser']
    user = config_dict['user']
    project = config_dict['project']
    run_id = config_dict['run_id']

    print("*****running!!!*****")
    
    srcname = cmd_line[1]
    Jfactor = 10**float(cmd_line[2]) #J-factor in units of GeV^2/cm^5
    sigmaJ = float(cmd_line[3]) #log10J_error


    sed_suffix = 'free_BG_sed'

    sedfits = '/zfs/astrohe/{}/Stacking_Analysis/{}/Preprocessed_Sources/{}/JLA/{}/output/{}_{}.fits'\
        .format(user, project, run_id,srcname, srcname, sed_suffix)

    save_array_path = '/zfs/astrohe/{}/Stacking_Analysis/{}/Numpy_Arrays/diMauro/{}/'\
                        .format(user, project, run_id)

    if os.path.exists(save_array_path) !=True:
        print('does not exist')
        print('creating path')
        os.system('mkdir -p '+save_array_path)
    else:
        print('path exists')


    # Inputs
    channelname = 'bb'
    EWcorr = 'No'
    indexDMchannel = 7


    print("Running for {}, J: {:.2e}, sigmaJ: {}".format(srcname,Jfactor, sigmaJ))


    sed = fits.open(sedfits)[1].data
    table_refdnde = copy.deepcopy(sed['ref_dnde'])
    table_normscan = copy.deepcopy(sed['norm_scan'])
    #convert table_normscan to units of cm^-2 s^-1 MeV^-1
    for t in range(len(table_refdnde)):
        table_normscan[t] = sed['norm_scan'][t,:]*table_refdnde[t]
        
    table_loglikescan = copy.deepcopy(sed['dloglike_scan'])

    # Convert MeV to GeV
    print("Define eref_vec in GeV...")
    eref_vec = copy.deepcopy(sed['e_ref'])/1.e3

    #-------------------------------------------
    # Producing Loglike profiling with sigmav.
    #-------------------------------------------

    # create results file
    print("create_vectors")
    mass_vec =  np.logspace(0,4,40)
    sigmav_vec = np.logspace(-28,-22,60)#np.logspace(-28,-23,50)
    LogLike_vec = np.zeros(shape=(len(mass_vec),len(sigmav_vec)))
    LogLikeJprof_vec = np.zeros(shape=(len(mass_vec),len(sigmav_vec)))

    for u in range(len(mass_vec)):
        for t in range(len(sigmav_vec)):
            DMprop = [mass_vec[u],sigmav_vec[t],indexDMchannel, EWcorr]
            LogLike_vec[u,t] = funcLogLike(DMprop,table_normscan,table_loglikescan, eref_vec=eref_vec, Jfactor=Jfactor)
            LogLikeJprof_vec[u,t] = funcLogLike_profJ(DMprop,table_normscan,table_loglikescan, eref_vec=eref_vec, Jfactor=Jfactor, sigmaJ=sigmaJ)[0]
            print('{:.2e} {:.2f}, {}'.format(sigmav_vec[t], mass_vec[u], np.round(LogLikeJprof_vec[u],2)))


    np.save(save_array_path+'{}_Jprior_freeBG_dlike.npy'.format(srcname), -LogLikeJprof_vec)
    np.save(save_array_path+'{}_noprior_freeBG_dlike.npy'.format(srcname), -LogLike_vec)


################################################
if __name__=="__main__":
    main(sys.argv)

