import pandas as pd
import matplotlib.pyplot as plt
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
import math
from os.path import exists

########################################################################################################################################################################
# Helper functions
########################################################################################################################################################################

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
    

def saveFromInput(saveFile, inputO:bool=False): 
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
		egDF.to_csv(saveFile, index=False)
	
	
########################################################################################################################################################################
# Overlap functions
########################################################################################################################################################################

def galacticPlaneOverlap(egDF, deg:float=15.):
#Remove candidates too close to galactic plane

	print("Removing EGs with b<15deg")
	removed=0
	toRemove=[]
	for x in egDF.index:
		if abs(egDF.loc[x,'b'])<deg:
			toRemove.append(x) #remove row if b<15deg
			removed+=1
	print(f"    Removed {removed} EGs - left with {(len(egDF.index)-removed)/len(egDF.index)*100:.2f}% of sample")
	return toRemove
	
def compCat(compDFName, raIn, decIn, egDF, sep:float=1.0):

	compDF=pd.read_csv(compDFName)
	toRemove=[]
	cl1=SkyCoord(l = egDF['l']*u.deg, b = egDF['b']*u.deg, frame='galactic').icrs
	cl2 = SkyCoord(ra = compDF[raIn]*u.degree,dec = compDF[decIn]*u.degree)

	#return arrays of elements from cl2 that match cl1
	#index of cl2 object closest to the corresponding cl1 object, 2D distance between objects, 3D distance between object (if exists)
	idx, d2d, d3d = cl1.match_to_catalog_sky(cl2)

	#Max separation between EG and other source
	max_sep = sep*u.degree #Maximum distance to neighbor
	sep_constraint = d2d < max_sep #array of booleans corresponding to cl1 - True if sep too small and cl1 element should be removed
	print(f"    Eliminating {sep_constraint.sum()} EGs for overlap within {sep} deg")

	for i in range(len(egDF)):
		if sep_constraint[i]:
			toRemove.append(i)
		
	return toRemove
	
########################################################################################################################################################################
# Main
########################################################################################################################################################################
	  
########################################################################################################################################################################
#Catalog from Kormendy and Ho (https://arxiv.org/pdf/1304.7762.pdf) with M_BH
print("Importing CSV 1 - Kormendy and Ho")
egDF=pd.read_csv('EGs_KH.csv')

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

#Remove candidates too close to galactic plane
toRemove = []
toRemove.append(galacticPlaneOverlap(egDF))

#Compare to BZCAT Blazar catalog (https://www.ssdc.asi.it/bzcat5/) 
print("Consider blazar catalog from BZCAT")
toRemove.append(compCat('bzcat_blazarCatalog.csv',' R.A. (J2000) ', ' Dec. (J2000) ', egDF))

#Compare to 2MRS radio galaxy catalog (http://ragolu.science.ru.nl/index.html) 
print("Consider radio galaxy catalog from 2MRS")
toRemove.append(compCat('2mrs_radioCatalog.csv','ra', 'dec', egDF))


#Remove EGs marked to remove
toRemove=sum(toRemove, []) #flatten list
toRemove = list(set(toRemove)) #remove duplicates
toRemove.sort()
for i in toRemove:
	egDF.drop(i,inplace=True)
egDF.reset_index(inplace=True)

print(f"Left with {len(egDF)} of the original {egDForigLength} EGs - {100*len(egDF)/egDForigLength:.2f}%")

#Save resulting EG list as a new csv
finalList="EGs_KH_overlapRemoved.csv"
saveFromInput(finalList,exists(finalList))
	
#Catalog from Kormendy and Ho (https://arxiv.org/pdf/1304.7762.pdf) with M_BH
########################################################################################################################################################################
