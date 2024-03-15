import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import glob,os,sys

#to run - send directory of csv's with spectrum info from "interpretation" code. Include ending /
	
############################################################################################################
#MAIN
############################################################################################################
def main():

	#Import CSVs
	df_spec=pd.DataFrame()
	df_spec_half=pd.DataFrame()
	df_spec_double=pd.DataFrame()
	
	for i,f in enumerate(glob.glob(sys.argv[1]+'/*.csv')):
		gal_name = f[len(sys.argv[1])+1:-4]
		df = pd.read_csv(f, names=["energy[GeV]", "flux"])
		if i==0:
			energy = pd.Series(df['energy[GeV]'])
		if 'nominal' in gal_name:
			df_spec[gal_name] = df['flux']
		elif 'half' in gal_name:
			df_spec_half[gal_name] = df['flux']
		else:
			df_spec_double[gal_name] = df['flux']
			
	"""	
	#Plot on one axis
	plt.plot(energy, df_spec, '-o')
	plt.legend(list(df_spec.columns), loc="best")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux [ph/cm^2/GeV]')
	plt.yscale('log')
	plt.show()
	
	plt.plot(energy, df_spec_half)
	plt.legend(list(df_spec_half.columns), loc="best")
	plt.title("Half distance")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux [ph/cm^2/GeV]')
	plt.yscale('log')
	plt.show()
	
	plt.plot(energy, df_spec_double)
	plt.legend(list(df_spec_double.columns), loc="best")
	plt.title("Double distance")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux [ph/cm^2/GeV]')
	plt.yscale('log')
	plt.show()
	"""
	
	#Compare individual galaxies 
	#4889 arbitrarily
	plt.plot(energy, df_spec['4889_nominal'], '-o', label = 'center')
	plt.plot(energy, [4.61900753e-28, 6.13759002e-27, 6.35972112e-26, 6.57952834e-25,
 6.77592157e-24, 6.87352663e-23, 6.65053099e-22, 5.57342377e-21,
 3.20343655e-20, 8.74126961e-20, 4.29060131e-20, 5.82836965e-21,
 3.97336498e-25], '-o', label='integrated')
	plt.legend(loc='best')
	plt.title("NGC4889")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux [ph/cm^2/GeV]')
	plt.yscale('log')
	plt.show()
	
	"""
	plt.plot(energy, df_spec['4889_nominal'], label='nominal')
	plt.plot(energy, df_spec_half['4889_half'], label='half')
	plt.plot(energy, df_spec_double['4889_double'], label='double')
	plt.legend(loc='best')
	plt.title("NGC4889")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux [ph/cm^2/GeV]')
	plt.yscale('log')
	plt.show()
	
	plt.plot(energy, df_spec['4889_nominal']/df_spec_double['4889_double'])
	plt.title("NGC4889 - nominal/double")
	plt.xlabel('Photon energy [GeV]')
	plt.xscale('log')
	plt.ylabel('DM flux ratio')
	plt.show()	
	"""

################################################################################################
################################################################################################
if __name__=="__main__":
	
	main()