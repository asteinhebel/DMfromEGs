import pandas as pd
import matplotlib.pyplot as plt
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
import math
from astropy.io import fits
from astropy.table import Table

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
    
########################################################################################################################################################################
#Catalog from Kormendy and Ho (https://arxiv.org/pdf/1304.7762.pdf) with M_BH
print("Importing CSV files")
egDF=pd.read_csv('EGs_new.csv')

#get l,b from Michela
egDF_mn=pd.read_csv('EGsTable.csv')
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

#Remove candidates too close to galactic plane
print("Removing EGs with b<15deg")
removed=0
toRemove=[]
for x in egDF.index:
	if abs(egDF.loc[x,'b'])<15:
		toRemove.append(x) #remove row if b<15deg
		removed+=1
print(f"    Removed {removed} EGs - left with {(len(egDF.index)-removed)/len(egDF.index)*100}% of sample")


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

#Compare to BZCAT Blazar catalog (https://www.ssdc.asi.it/bzcat5/) 
print("Consider blazar catalog from BZCAT")
blazarDF=pd.read_csv('bzcat_blazarCatalog.csv')

cl1 = c.icrs #egDF SkyCoord object
cl2 = SkyCoord(ra = blazarDF[' R.A. (J2000) ']*u.degree,dec = blazarDF[' Dec. (J2000) ']*u.degree)

print ("Comparing EGs to BZCAT")

#return arrays of elements from cl2 that match cl1
#index of cl2 object closest to the corresponding cl1 object, 2D distance between objects, 3D distance between object (if exists)
idx, d2d, d3d = cl1.match_to_catalog_sky(cl2)

#Max separation between EG and other source
sep=1.0
max_sep = sep*u.degree #Maximum distance to neighbor
sep_constraint = d2d < max_sep #array of booleans corresponding to cl1 - True if sep too small and cl1 element should be removed
print(f"    Eliminating {sep_constraint.sum()} EGs for overlap with BZCAT catalog - within {sep} deg")

for i in range(len(egDF)):
	if sep_constraint[i]:
		toRemove.append(i)

#Compare to 2MRS radio galaxy catalog (http://ragolu.science.ru.nl/index.html) 
print("Consider radio galaxy catalog from 2MRS")
radioDF=pd.read_csv('2mrs_radioCatalog.csv')
cl3 = SkyCoord(ra = radioDF['ra']*u.degree,dec = radioDF['dec']*u.degree)
print ("Comparing EGs to 2MRS")
idx, d2d, d3d = cl1.match_to_catalog_sky(cl3)
max_sep = sep*u.degree #Maximum distance to neighbor
sep_constraint = d2d < max_sep #array of booleans corresponding to cl1 - True if sep too small and cl1 element should be removed
print(f"    Eliminating {sep_constraint.sum()} EGs for overlap with 2MRS catalog - within {sep} deg")

for i in range(len(egDF)):
	if sep_constraint[i]:
		toRemove.append(i)


#Remove EGs marked to remove
toRemove=list(set(toRemove))
toRemove.sort()
for i in toRemove:
	egDF.drop(i,inplace=True)

print(f"Left with {len(egDF)} of the original {egDForigLength} EGs - {100*len(egDF)/egDForigLength}%")

#Save resulting EG list as a new csv
finalList="EGs_overlapRemoved.csv"
print(f"Saving the remaining EGs to {finalList}")
egDF.to_csv(finalList, index=False)
