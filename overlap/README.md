# Identify viable elliptical galaxy (EG) sources for Fermi analysis

Suitable sources for analysis must:
- Be >10deg above the galactic plane
- Not fall within 0.1deg of a known radio galaxy (from 2MRS catalog)
- Not fall within 0.1deg of a known blazar (from BZCAT catalog)
- Not fall within 0.1deg of a known gamma ray source (from Fermi 4FGL-2DL catalog)

## Two input EG catalogs:
- ``KH" from Kormendy and Ho (<https://arxiv.org/pdf/1304.7762.pdf>)
- ``DB" from Dabringhausen and Fellhauer (<https://academic.oup.com/mnras/article/460/4/4492/2609151>)

Initial EG catalogs must be saved as CSV and contain at least the variables:
1. Index number (index)
2. NGC name (Name)
3. Distance (Distance [Mpc])
4. RA (RA)
5. Dec (Dec)
6. Galactic longitude (l)
7. Galactic latitude (b)
8. Stellar velocity dispersion (sigma_e [km/s])

For KH, this is done in catalogOverlap.py code by adding values from the supplementary EGs_KH_lb.csv 
For DF, this is done in the creation of the input EGs_DF.csv. To create this original input from the article's Tables (saved as .txt files).

	cd tables_DF
	
	python makeCSV.py 
Will output EGs_DF.csv that can be copied to the directory above.

## To run overlap code 

### In terminal:
	python catalogOverlap.py
	
### As iPython notebook:
	jupyter notebook
And open catalogOverlap.ipynb

### Details
Code requires that original csv's and csvs from overlap catalogs be stored in the same directory. This includes:
- EGs_DF.csv
- EGs_KH.csv
- EGs_KH_lb.csv
- 2mrs_radioCatalog.csv
- bzcat_blazarCatalog.csv
- fermi_4fgl_gammaCatalog.csv

Code outputs:
- EGs_DF_overlapRemoved.csv
- EGs_KH_overlapRemoved.csv

If one of these outputs already exists, code will prompt in-terminal whether the original file should be overwritten or not.

Runtime optional arguments:
- -o, directory for saving outputs, needs string input
- -p, plot histograms of minimum separation between EG and comparison catalog entries
- -s, save output csv lists of good EGs and separation histograms 