import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy.table import Table
import math
import os
import argparse

########################################################################################################################################################################
# Helper functions
########################################################################################################################################################################

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
    

def saveFromInput(saveFile, df, inputO:bool=False): 
	if inputO:
		inp= 'a'
		while inp!='y' and inp!='n':
			print("Overwrite? y/n")
			inp = input()
		savePlt=True if inp=='y' else False	
	else:
		savePlt=True
	
	if savePlt: 
		print(f"Saving the remaining EGs to {saveFile}")
		df.to_csv(saveFile, index=False)

def addMBH():
	pass	
	
########################################################################################################################################################################
# Overlap functions
########################################################################################################################################################################
def extraKHprocessing(egDF):
	#get l,b from Michela
	egDF_mn=pd.read_csv('EGs_KH_lb.csv')
	egDF_mn[['Name','Type','Distance[Mpc]','Distance Uncert','Ks','M_KsT','M_VT','(V-Ks)_0','(B-V)_0','M_BH[Msolar*10e6]','Flags:M','Flags:C','Flags:M_BH','l','b']]=egDF_mn['Name;Type;Distance[Mpc];Ks;M_KsT;M_VT;(V-Ks)_0;(B-V)_0;M_BH[Msolar*10e6];Flags:M;Flags:C;Flags:M_BH;;l;b'].str.split(';', expand=True)
	egDF_mn.drop(columns=['Name;Type;Distance[Mpc];Ks;M_KsT;M_VT;(V-Ks)_0;(B-V)_0;M_BH[Msolar*10e6];Flags:M;Flags:C;Flags:M_BH;;l;b'], inplace=True)
	egDF['l']=egDF_mn['l']
	egDF['b']=egDF_mn['b']
	
	#Remove rows of NaN entries
	egDF.dropna(inplace=True)
	egDForigLength=len(egDF)

	#Convert l,b to floats
	toConvert=['l','b']
	for colu in toConvert:
		egDF[colu]=pd.to_numeric(egDF[colu])

	#Add RA/Dec in degrees to egDF
	c=SkyCoord(l = egDF['l']*u.deg, b = egDF['b']*u.deg, frame='galactic')

	eg_ra=[]
	eg_dec=[]
	for coord in c.icrs: #convert galactic coordinates to icrs RA/Dec
		strcoord=coord.to_string()
		eg_ra.append(float(strcoord.split(' ')[0]))
		eg_dec.append(float(strcoord.split(' ')[1]))
	egDF['RA']=eg_ra
	egDF['Dec']=eg_dec
	
	return egDF

def galacticPlaneOverlap(egDF, deg:float=10.):
	"""Remove candidates too close to galactic plane
	Return: Array of index locations of EGs that should be removed due to galactic plane proximity"""

	print(f"Removing EGs with b<{deg} deg")
	removed=0
	toRemove=[]
	for x in egDF.index:
		if abs(egDF.loc[x,'b'])<deg:
			toRemove.append(x) #remove row if b<15deg
			removed+=1
	print(f"    Removed {removed} EGs - left with {(len(egDF.index)-removed)/len(egDF.index)*100:.2f}% of sample")
	return toRemove
	
def compCat(compDFName, raIn, decIn, egDF, egDFName, sep:float=0.1):
	"""Identify EG candidates that are too close to sources in another catalog
	Return: Array of index locations of EGs that should be removed due to proximity to item in other catalog"""

	compDF=pd.read_csv(compDFName)
	toRemove=[]
	cl1=SkyCoord(l = egDF['l']*u.deg, b = egDF['b']*u.deg, frame='galactic').icrs
	cl2 = SkyCoord(ra = compDF[raIn]*u.degree,dec = compDF[decIn]*u.degree)

	#return arrays of elements from cl2 that match cl1
	#index of cl2 object closest to the corresponding cl1 object, 2D distance between objects, 3D distance between object (if exists)
	idx, d2d, d3d = cl1.match_to_catalog_sky(cl2)
	
	if args.plot:
		#Make histograms of minimum distance between each EG and the closest entry from the comparison catalog
		degMax = np.max(d2d.degree)
		binEdges = np.arange(0, degMax, 0.05) #0.05deg bins
		np.concatenate((binEdges,[degMax]))
		hist=plt.hist(d2d.degree, bins=binEdges)
		plt.title(compDFName[:-4])
		plt.xlabel("Min separation distance [deg]")
		plt.ylabel("elliptical galacy count")
		if args.saveOut:
			savename=f"{args.outDir}/hist_minSep_{egDFName}_{compDFName[:-4]}.png"
			print(f"Saving {savename}")
			plt.savefig(savename)
		else:
			plt.show()
		plt.clf()

	#Max separation between EG and other source
	max_sep = sep*u.degree #Maximum distance to neighbor
	sep_constraint = d2d < max_sep #array of booleans corresponding to cl1 - True if sep too small and cl1 element should be removed
	print(f"    Eliminating {sep_constraint.sum()} EGs for overlap within {sep} deg")

	print("EGs TO BE REMOVED:")
	to_remove = egDF['Name']*(d2d<max_sep)
	to_remove.replace('', np.nan, inplace=True)
	to_remove.dropna(inplace=True)
	print(to_remove)

	for i in range(len(egDF)):
		if sep_constraint[i]:
			toRemove.append(i)
		
	return toRemove
	
########################################################################################################################################################################
# Main
########################################################################################################################################################################
	  
def main(args, f_in, extraprocess): 

	egDF=pd.read_csv(f_in)
	f_in_name=f_in[:-4]
	
	if extraprocess:
		egDF=extraKHprocessing(egDF)

	egDForigLength=len(egDF)

	#Remove candidates too close to galactic plane
	toRemove = []
	toRemove.append(galacticPlaneOverlap(egDF))

	#Compare to BZCAT Blazar catalog (https://www.ssdc.asi.it/bzcat5/) 
	print("Consider blazar catalog from BZCAT")
	toRemove.append(compCat('bzcat_blazarCatalog.csv',' R.A. (J2000) ', ' Dec. (J2000) ', egDF, f_in_name))

	#Compare to 2MRS radio galaxy catalog (http://ragolu.science.ru.nl/index.html) 
	print("Consider radio galaxy catalog from 2MRS")
	toRemove.append(compCat('2mrs_radioCatalog.csv','ra', 'dec', egDF, f_in_name))

	#Compare to 4FGL gamma-ray source catalog (https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/) 
	print("Consider 4FGL gamma-ray catalog")
	toRemove.append(compCat('fermi_4fgl_gammaCatalog.csv','RAJ2000', 'DEJ2000', egDF, f_in_name))

	#Remove EGs marked to remove
	toRemove=sum(toRemove, []) #flatten list
	toRemove = list(set(toRemove)) #remove duplicates
	toRemove.sort()
	for i in toRemove:
		egDF.drop(i,inplace=True)
	egDF.reset_index(inplace=True)

	print(f"Left with {len(egDF)} of the original {egDForigLength} EGs - {100*len(egDF)/egDForigLength:.2f}%")

	if args.saveOut:
		#Save resulting EG list as a new csv
		finalList=f"{args.outDir}/{f_in_name}_overlapRemoved.csv"
		saveFromInput(finalList,egDF,os.path.exists(args.outDir+"/"+finalList))		


########################################################################################################################################################################
# call main
########################################################################################################################################################################
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Determine overlap of source EGs with other catalogs')
	parser.add_argument('-s', '--saveOut', action='store_true', default=False, required=False, 
		help='Save all plots (no display_ and/or csv lists to current dir. Default: False (only displays/prints to terminal, does not save)')
	parser.add_argument('-p', '--plot', action='store_true', default=False, required=False, 
		help='Plot histograms of distance between EG and closest entry in other catalogs. Default:False')
	parser.add_argument('-o', '--outDir', default='./', required=False,  
        help='Directory to save outputs. Default: current directory')

	parser.add_argument
	args = parser.parse_args()
	
	#create output dir if does not exist
	if not os.path.exists(args.outDir):
		os.mkdir(args.outDir)
	
	extraprocessing = [True, False]
	f_in = ["EGs_KH.csv","EGs_DF.csv"]
	for i,extra in enumerate(extraprocessing):
		print(f"Importing CSV {i+1} - {f_in[i]}")
		main(args, f_in[i], extra)
			
	#Catalog KH from Kormendy and Ho (https://arxiv.org/pdf/1304.7762.pdf) with M_BH
	#Catalog DF from Dabringhausen and Fellhauer (https://academic.oup.com/mnras/article/460/4/4492/2609151) without M_BH	