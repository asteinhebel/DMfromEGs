# -*- coding: utf-8 -*-
#imports
import os,sys
from math import *
import numpy as np
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import copy
import astropy.io.fits as fits
from scipy import interpolate



def plotTS(summed, masses, sigmav, vmin=0, interp=True,dof=2,
           save=False, filename='TS.png', title=None):

    fig = plt.figure(figsize=(3.5,3),dpi=150)
    fig.gca().patch.set_facecolor('white')
    plt.loglog()
    
    if interp:
        f = interpolate.RectBivariateSpline(sigmav,masses,summed)

        mass = np.logspace(np.log10(masses.min()),np.log10(masses.max()),num=200,endpoint=True)
        cross = np.logspace(np.log10(sigmav.min()),np.log10(sigmav.max()),num=200,endpoint=True)
        ts = f(cross, mass)
        img = plt.pcolormesh(mass,cross,ts,cmap="inferno",vmin=vmin,vmax=summed.max())
    else:
        img = plt.pcolormesh(masses,sigmav,summed,cmap="inferno",vmin=vmin,vmax=summed.max())

    ind = np.unravel_index(np.argmax(summed,axis=None),summed.shape)

    ind = np.unravel_index(np.argmax(summed,axis=None),summed.shape)


    best_index_value = ind[0]
    best_flux_value = ind[1]

    best_index = sigmav[ind[0]]
    best_flux = masses[ind[1]]

    if(dof!=None):
        if dof==2:
            levels = summed.max() - np.asarray([11.8,6.17,2.3])
        if dof==3:
            levels = summed.max() - np.asarray([14.16, 8.02, 3.53])
        plt.contour(masses,sigmav,summed,levels=levels,colors='limegreen',linestyles=["-.",'--',"-"], linewidths=3*[1],alpha=1)

    plt.plot(best_flux, best_index,marker="+",ms=4,color="black")
    
    
    if(title!=None):
    	plt.title(title, fontsize=10)
    plt.xlabel('E [GeV]')
    plt.ylabel('E$^2$dN/dE [MeV/cm$^2$/s]')
    cbr = plt.colorbar(img, shrink=1)
    cbr.ax.set_title('TS')
    plt.tight_layout()
    
    
    print(f"Peak TS: {summed.max():.2f}")
    print(f"Peak energy: {best_flux:.3f} GeV")
    print(f"Peak scan bin : {best_index} MeV/cm$^2$/s")
    
    
    if(save):
        plt.savefig(filename)
    else:
    	plt.show()
    	
    return best_index, best_flux, summed.max()

############################################################################################################
#MAIN
############################################################################################################
def main(cmd_line):

	srcname = sys.argv[1]
	sedfits =  f"{homepath}{srcname}/output/{srcname}_free_BG_sed.fits"
	save_array_path = f"{homepath}{srcname}/output/plots/"
	sed = fits.open(sedfits)[1].data
	like_file = copy.deepcopy(sed['dloglike_scan'])
	
	en_vec = copy.deepcopy(sed['e_ref'])/1.e3
	flux_vec = np.logspace(-8,-5,20)
	
	#correct 0 point
	TS_array = 2*(like_file-like_file[0,0])
	
	plotTS(TS_array.T, en_vec, flux_vec, interp=False, vmin=0, save=savePlots, filename=homepath+srcname+"/output/plots/TS_enFlux.png", title=srcname)


################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":
	homepath = '/Users/asteinhe/FermiLAT/BHinEGs_DM/run/'
	savePlots=False
	main(sys.argv)