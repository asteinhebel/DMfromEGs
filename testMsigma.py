import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
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

#Parse sigma and errors from CSV string
def getSigma(sigmaIn):
	sig, errUp, errDown = [], [], []
	for sigma in sigmaIn:
		try:
			dataIn = sigma.split('±')
			sig.append(float(dataIn[0]))
			errUp.append(float(dataIn[1]))
			errDown.append(float(dataIn[1]))
		except ValueError: #unequal uncertainties for plus/minus
			dataIn = sigma.split('+')
			if dataIn[0]!= '...':
				sig.append(float(dataIn[0]))
				errIn = dataIn[1].split('-')
				errUp.append(float(errIn[0]))
				errDown.append(float(errIn[1]))
			else:
				sig.append(np.NaN)
				errUp.append(np.NaN)
				errDown.append(np.NaN)
	return sig, errUp, errDown

#Use M-sigma relation to calculate M from sigma
def calcMBH(sigma, a:float=-0.501, b:float=4.377):
	#Defaults from KH, Sec. 6.6, Eqn. 3
	m = math.exp(a) * (sigma/200.)**b * 1000.
	return m	
	
#Parse measured BH mass errors from CSV string
def getBHErrs(massIn):
	errUp, errDown = [], []
	for mass in massIn:
		mass = mass.replace('(',';').replace(')',';')
		dataIn = mass.split(';')
		meas= float(dataIn[0])
		errs = dataIn[1]
		scale = float(dataIn[2][-1])
		if scale>1:
			exp = scale-6
		else: #10^10, overflowed single digit scale counter
			exp = 4. #(10^10 scaled by 10^6)
		#CSV gives low/high so actual difference value from measured must be calculated
		errUp.append((float(errs.split('−')[1])-meas)*10**exp)
		errDown.append((meas-float(errs.split('−')[0]))*10**exp)
	return errUp, errDown
	
	
########################################################################################################################################################################
# Main
########################################################################################################################################################################
	  
### ATTENTION USER!!! ###
# Choose which csv file to run over
egDF=pd.read_csv('EGs_KH.csv')
#egDF=pd.read_csv('EGs_KH_overlapRemoved.csv')

#Define directory to save plots to
outDir = 'plots_testMsigma/'

#Remove csv entries with NaNs
egDF.dropna(inplace=True)

#Calculate MBH from sigma
calc1, calc2 = [], []
meas = egDF['M_BH [Msolar * 10e6]']
ind = np.arange(len(meas))
massErrUp, massErrDown = getBHErrs(egDF['M_BH (low - high) [Msolar]'])
sig, sigErrUp, sigErrDown = getSigma(egDF['sigma_e [km/s]'])
for i,sigma in enumerate(sig):
	if sigma!='...':
		calc1.append(calcMBH(sigma)) #Defaults from KH, (full sample fit) Sec. 6.6, Eqn. 3 
		calc2.append(calcMBH(sigma,-0.54,4.26))#KH, Sec. 6.7.2 (low-sigma fit only), Eqn. 12
	else: #M_BH does not exist in orig dataset
		calc1.append(0.)
		calc2.append(0.)

#Find necessary error to include on calculated points
np_meas = np.array(meas)
np_calc = np.array(calc1)
np_relDif = np.abs((np_calc-np_meas)/np_calc)
plt.hist(np_relDif)
plt.xlabel('Absolute relative difference between M-sigma mass and measured')
plt.ylabel('Counts')
savePlot = 'hist_calcErr'
plt.savefig(f"{outDir}{savePlot}.png")
plt.show()

#Define error to introduce to calculated M_BH - attribute to originating EG
maxIndex = np.argmax(np_relDif)+1
maxVal = max(np_relDif)
print(f"Max relative difference = {maxVal.round(2)}")
print(f"Max value at index {maxIndex} - EG {egDF['Name'].iloc[maxIndex]}")
print(f"sigma {egDF['sigma_e [km/s]'].iloc[maxIndex]}")
print(f"M_BH {egDF['M_BH (low - high) [Msolar]'].iloc[maxIndex]}")
print(f"Measured M_BH {meas[maxIndex]}")
print(f"Calculated M_BH {calc1[maxIndex]}")

#Define error arrays
sigErr = [sigErrUp, sigErrDown]
massErr = [massErrUp, massErrDown]
calcErr=np_calc*maxVal

#Plot M_BH vs sigma with uncertainty on sigma and only measured M_BH - calculated M_BH and measured. Vertical bar indicates sigma of EG with the largest disagreement between measured and calculated M_BH
plt.vlines(x=sig[maxIndex], ymin=1, ymax=100000, colors='red')
plt.errorbar(sig, calc1, xerr=sigErr, label='Calculated', fmt='.', color='b')
plt.errorbar(sig, calc2, xerr=sigErr, label='Calculated, low sig', fmt='.', color='y')
plt.errorbar(sig, meas, xerr=sigErr, yerr=massErr, label='Measured', fmt='.', color='g')
plt.xlabel('sigma (stellar velocity dispersion) [km/s]')
plt.ylabel('M_BH * 10^6 M_solar')
plt.title('Breakdown of M-sigma relation when sigma>270 km/s')
plt.legend(loc='best')
plt.yscale('log')
savePlot = 'sigmaVSmbh_noCalcErr'
plt.savefig(f"{outDir}{savePlot}.png")
plt.show()

#Plot M_BH vs sigma - same as plot above but with relative errors added to calculated M_BH 
plt.errorbar(sig, calc1, xerr=sigErr, yerr=calcErr, label='Calculated', fmt='.', color='b')
plt.errorbar(sig, calc2, xerr=sigErr, label='Calculated, low sig', fmt='.', color='y')
plt.errorbar(sig, meas, xerr=sigErr, yerr=massErr, label='Measured', fmt='.', color='g')
plt.xlabel('sigma (stellar velocity dispersion) [km/s]')
plt.ylabel('M_BH * 10^6 M_solar')
plt.title('Breakdown of M-sigma relation when sigma>270 km/s')
plt.legend(loc='best')
plt.yscale('log')
savePlot = 'sigmaVSmbh_calcErr'
plt.savefig(f"{outDir}{savePlot}.png")
plt.show()

#Plot sigma vs relative difference between measured and calculated M_BH to see if calculated M_BH differs as a function of sigma
plt.errorbar(sig, np_relDif, xerr=sigErr, fmt='o')
plt.xlabel('sigma (stellar velocity dispersion) [km/s]')
plt.ylabel('Absolute relative difference between M-sigma mass and measured')
savePlot = 'sigmaVScalcErr'
plt.savefig(f"{outDir}{savePlot}.png")
plt.show()