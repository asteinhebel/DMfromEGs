import pandas as pd
import numpy as np
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from os.path import exists

# Necessary Columns:
## Index number (index)
## NGC name (Name)
## Distance (Distance [Mpc])
## RA (RA)
## Dec (Dec)
## Galactic longitude (l)
## Galactic latitude (b)
## Stellar velocity dispersion (sigma_e [km/s])

###########################################################################################
# Helper functions
###########################################################################################

def getCoords(df):

	ra = df2['asc']
	ra_str = []
	for val in ra:
		c1 = val.replace(":","h",1)
		c2 = c1.replace(":","m",1)
		c3 = c2+"s"
		ra_str.append(c3)
		
	dec = df2['dec']
	dec_str = []
	for val in dec:
		d1 = val.replace(":","d",1)
		d2 = d1.replace(":","m",1)
		d3 = d2+"s"
		dec_str.append(d3)

	ra_np = np.array(ra_str)
	dec_np = np.array(dec_str)
	pairs = np.apply_along_axis(' '.join, 0, [ra_np, dec_np])

	return pairs
	
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
		print(f"Saving to {saveFile}")
		df.to_csv(saveFile, index=False)

###########################################################################################
# DF1 - get index, name
###########################################################################################
df1 = pd.read_fwf('Table1(names).txt', sep=" ", skiprows=25, dtype=str, header=None)
cols1 = ["index", "NGC", "IC", "UGC","AGC","PGC","CGCG ","VCC","GMP","FCC","LLH2012","MHM2009","HCC","Abell_496", "Other"]
df1.columns = cols1

# If 'other' name exists, use that name. Then use NGC, then IC, etc from cols array
## In code, iterate backwards through cols array so each item is overridden. Put 'other' still at end of array so that always overrides if it exists
nm = [''] * len(df1)
labels=cols1[1:-1]
labels.reverse()
keysForNaming = [labels,[cols1[-1]]]
keysForNaming = sum(keysForNaming, [])

for i in range(len(df1)):
	for k in keysForNaming:
		if df1[k].iloc[i]!="0" and k!='Other':
			nm[i]=(k+" "+df1[k].iloc[i])
		elif df1[k].iloc[i]!="0" and k=='Other':
			nm[i]=(df1[k].iloc[i])
	#If no name in any category, entry remains empty (none should fall in this category)
		
#Add name and index to output DF		
dfOut = pd.DataFrame(nm,columns =['Name'])
ind=df1['index']
dfOut.insert(0, "index", ind)

###########################################################################################
# DF2 - get RA, Dec, distance
###########################################################################################
df2 = pd.read_fwf('Table2(locations).txt',skiprows=20, header=None, sep=" ", dtype=str)
cols2 = ["id","asc","dec","s.","d","N_d","s.","v_rad","s."]
df2.columns = cols2

#Add distance to output DF
distance = df2['d']
dfOut['Distance [Mpc]'] = pd.to_numeric(distance)

#Add RA, Dec, l, b to output DF
eg_ra, eg_dec = [], []
eg_l, eg_b = [], []

#Convert RA and Dec from J2000.0 h:m:s / d:m:s to deg
strCoords = getCoords(df2)
c=SkyCoord(strCoords)

for coord in c.icrs: 
	strcoord=coord.to_string()
	#save ra/dec
	eg_ra.append(float(strcoord.split(' ')[0]))
	eg_dec.append(float(strcoord.split(' ')[1]))
	#save l/b
	galcoord=coord.galactic.to_string()
	eg_l.append(float(galcoord.split(' ')[0]))
	eg_b.append(float(galcoord.split(' ')[1]))

dfOut['RA']=eg_ra
dfOut['Dec']=eg_dec
dfOut['l']=eg_l
dfOut['b']=eg_b

###########################################################################################
# DF3 - get velocity dispersion
###########################################################################################
df3 = pd.read_fwf('Table11(structural-properties).txt',skiprows=33, header=None, sep=" ")

sigma_e = df3[7]
dfOut['sigma_e [km/s]'] = sigma_e

###########################################################################################
# Save output CSV
###########################################################################################

csvOut="EGs_DF.csv"
saveFromInput(csvOut,dfOut,exists(csvOut))
