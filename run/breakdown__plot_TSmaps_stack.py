# -*- coding: utf-8 -*-

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import os, sys, shutil

############################################################################################################
# HELPER FUNCTIONS
############################################################################################################

def plotTS(summed, masses, sigmav,vmin=0, interp=True,dof=2,save=False, filename='TS.png', title=None):

	fig = plt.figure(figsize=(3.5,3),dpi=150)
	fig.gca().patch.set_facecolor('white')
	plt.loglog()
	if interp:

		f = interpolate.RectBivariateSpline(sigmav,masses,summed)

		mass = np.logspace(np.log10(masses.min()),np.log10(masses.max()),num=200,endpoint=True)
		cross = np.logspace(np.log10(sigmav.min()),np.log10(sigmav.max()),num=200,endpoint=True)
		ts = f(cross, mass)
		img = plt.pcolormesh(mass,cross,ts,cmap="inferno",vmin=vmin,vmax=summed.max()*1.2)
	else:
		img = plt.pcolormesh(masses,sigmav,summed,cmap="inferno",vmin=vmin,vmax=summed.max())

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
	
	"""
	#DEBUG - FOR TESTING - OVERLAY LAT DATA FROM 2108.13646 
	latDF = pd.read_csv("latData.csv")
	latDF.columns = ['x','y']
	plt.plot(latDF['x'], latDF['y'], '-', label = "LAT")
	"""
	
	if(title!=None):
		plt.title(title, fontsize=10)
	plt.xlabel('$M_{\chi}$ [GeV]')
	plt.ylabel(r'$\left<\sigma v\right>$ [cm$^{3}$ s$^{-1}$]')
	cbr = plt.colorbar(img, shrink=1)
	cbr.ax.set_title('TS')
	plt.tight_layout()
	
	print(f"Peak TS: {summed.max():.2f}")
	print(f"Peak mass: {best_flux:.3f} GeV")
	print(f"Peak <sigmav> : {best_index} cm^3/s")
	
	if(save):
		plt.savefig(filename)
		with open(f'{filename[:-4]}.txt', 'w') as f:
			f.write(f"Peak TS: {summed.max():.2f} \n")
			f.write(f"Peak mass: {best_flux:.3f} GeV \n")
			f.write(f"Peak <sigmav> : {best_index} cm^3/s")
		f.close()
	else:
		plt.show()
		
	return best_index, best_flux, summed.max()
	
def makeDir(dirpath):
	if(os.path.isdir(dirpath)==True):
		print(f"Path to {dirpath} exists")
	else:
		print(f"path to {dirpath} does not exist, creating...")
		os.system('mkdir -p %s' %dirpath)
	return dirpath

############################################################################################################
#MAIN
############################################################################################################
def main(cmd_line):

	dmMassRange = [int(cmd_line[2]), int(cmd_line[3]), int(cmd_line[4])]
	dmXSRange = [int(cmd_line[5]), int(cmd_line[6]), int(cmd_line[7])]
	mass_vec = np.logspace(dmMassRange[0], dmMassRange[1], dmMassRange[2])
	sigmav_vec = np.logspace(dmXSRange[0], dmXSRange[1], dmXSRange[2])
	summed=0.
	
	#get sources
	print(f"Reading from {cmd_line[1]}")
	srclist = open(cmd_line[1],'r').read().split('\n')
	#protect against empty lines
	if '' in srclist:
		srclist.remove('')
	
	for i, srcname in enumerate(srclist):
		# The filename in the load functions need to be your specific numpy arrays for the run
		array_path = homepath+srcname+'/output/dloglike/'+subdir
		try:
		  like_file = np.load(array_path+'/{}_dlike.npy'.format(srcname))
		  print(f"Opening {srcname}")
		except:
		  print("{} Does not exist".format())
		  
		#Define TS
		TS_array = -2*like_file
		
		"""
		np.savetxt("TEST.csv", TS_array, delimiter=",")
		plt.hist(TS_array)
		plt.show()
		plt.clf()
		"""
		
		#Add to summed plot
		summed+= TS_array
		
		#save TS plots for individual inputs
		plotTS(TS_array.T, mass_vec, sigmav_vec, vmin=0, save=savePlots, filename=makeDir(homepath+srcname+"/output/plots/"+subdir+"/")+"TS.png", title=srcname)
	
	#save TS plot for stacked sample
	plotTS(summed.T, mass_vec, sigmav_vec, vmin=0, save=savePlots, filename=makeDir(homepath+"stack/"+subdir+"/")+"TS.png", title="stack")

################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":

	homepath = '/Users/asteinhe/FermiLAT/BHinEGs_DM/run/'
	subdir = 'mbh1'
	savePlots = True
	
	main(sys.argv)