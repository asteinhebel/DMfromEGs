#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import sys, os, shutil
import astropy.io.fits as fits
import copy

############################################################################################################
# HELPER FUNCTIONS
############################################################################################################
def makeDir(dirpath):
	if(os.path.isdir(dirpath)==True):
		print(f"Path to {dirpath} exists")
	else:
		print(f"path to {dirpath} does not exist, creating...")
		os.system('mkdir -p %s' %dirpath)
	
def diffPhotSpectrum(m, sigmav, mBH, dist, eref): 
	#DM spectrum
	sig26 = sigmav * 1e26 #scale to get sigma ^-26
	#Calculate dE/dN of DM spectrum, units of ph/cm^2/s/GeV (after multiplication with DM sigma), return as array of points that span LAT sensitivity 
	norm = (2.e-20)*pow(mBH,(5./2.))*sig26*pow(m,-2.)*pow(dist,-2.) * 2. / m
	spectrum_dm = np.array([ norm * (en/m)**3 * np.exp(-(en/m)**2) for en in eref])
	interp_spectrum_dm = interpolate.InterpolatedUnivariateSpline(eref, spectrum_dm, k=1, ext=0) #linear interpolation because k=1, ext=0=>extrapolate. Identical to interp1d

	#MSP spectrum
	spectrum_msp = np.array([1.2e-9 * (mBH / 1000) * pow(dist,-2.) * pow(en, -1.29) * np.exp(-en/2.5) for en in eref])
	interp_spectrum_msp = interpolate.InterpolatedUnivariateSpline(eref, spectrum_msp, k=1, ext=0)
	
	#Add DM + MSP
	spectrum = spectrum_dm + spectrum_msp
	interp_spectrum = interpolate.InterpolatedUnivariateSpline(eref, spectrum, k=1, ext=0)
	

	"""
	#DEBUG
	if m==100 and sigmav == 1.3141473626117527e-23:
		plt.plot(eref, spectrum, "o", label="DM")
		plt.plot(eref, interp_spectrum(eref), "-", label="DM interp")
		plt.plot(eref, spectrum_msp, "o", label="MSP")
		plt.plot(eref, interp_spectrum_msp(eref), "-", label="MSP interp")
		plt.legend(loc='best')
		plt.xlabel("Center of SED bins [GeV]")
		plt.ylabel("Expected flux [ph/cm^2/GeV]")
		plt.xscale('log')
		plt.yscale('log')
		plt.show()
		plt.clf()
	
		#save spectrum as csv for later plotting
		spec_out = np.asarray([eref, spectrum])
		np.savetxt(save_array_path+'spectrum_'+str(m)+'_'+str(sigmav)+'.csv', spec_out.T, delimiter=",")
		spec_out_msp = np.asarray([eref, spectrum_msp])
		np.savetxt(save_array_path+'spectrum_msp_'+str(m)+'_'+str(sigmav)+'.csv', spec_out_msp.T, delimiter=",")
	"""
	
	#return spectrum
	return interp_spectrum

def flux_DM_model(DMprop, spec, emin, emax, eref):
	#Compute dNdE DM from model
	DMmass = DMprop[0]
	sigmav = DMprop[1]	
	mBH = DMprop[2]
	dist = DMprop[3]
	
	"""
	#DEBUG
	#plt.plot(eref, fullSpec, "o")
	plt.plot(eref, spec(eref), "-")
	plt.xlabel("Observed photon energy [GeV]")
	plt.ylabel("Expected DM flux [ph/cm^2/s/GeV]")
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	plt.clf()
	"""
	
		
	#Integrate photon flux in each SED bin from min-max value (emin to emax)
	fluxInBins=[]
	for sedbin in range(len(emin)):
		#define finer mesh of points within SED bin to compute integral
		xx = np.linspace(emin[sedbin], emax[sedbin], 200)
		#integrate interpolated flux over this range
		fluxInBins.append(spec.integral(xx[0], xx[-1]))
		"""
		if DMprop[0]==100 and DMprop[1] == 1.3141473626117527e-23:
			print(f"Bin {sedbin}: [{emin[sedbin]}, {emax[sedbin]}]") 
			print(f"Bin center flux value = {spec(eref[sedbin])}")
			print(f"Bin edge flux values: {spec(emin[sedbin])}, {spec(emax[sedbin])}")
			print(f"Integrate from {xx[0]} to {xx[-1]}")
			print(f"Integral = {fluxInBins[-1]}")
		"""

	fluxInBins = np.array(fluxInBins)
	
	return fluxInBins
	
	
def funcLogLike(DMprop, spec, table_norm, table_loglike, emin, emax, eref):
	
	LogLike = []
	fluxval = flux_DM_model(DMprop, spec, emin, emax, eref) 
	"""
	if DMprop[0]==100 and DMprop[1] == 1.3141473626117527e-23:
		print("integrated flux vals")
		print(fluxval)
	"""
	#array the same length as number of SED bins, each value has units of ph/cm^2/s		
	
	#Do not consider first point of each scan - normalization of zero 
	table_norm = table_norm[:, 1:]
	
	table_loglike = table_loglike[:, 1:]
		
	#print(f"IN FUNCLOGLIKE WITH DM MASS {DMprop[0]} AND DM SIGMAV {DMprop[1]}")
	#calculate dloglike of DM flux value individually in each energy bin from interpolated data
	for t in range(len(table_norm)):
		#determine whether bin is suitable for interpolating (if dloglike max at lowest flux value)
		maxBin = np.argmax(table_loglike[t])
		if maxBin==0 and fluxval[t] < max(table_norm[t]):
		#if fluxval[t]>0. and fluxval[t] < max(table_norm[t]):
			#approximate function SED loglike = y = f(x) = f(table_norm)
			loglike_interp = interpolate.InterpolatedUnivariateSpline(table_norm[t], table_loglike[t], k=1, ext=0) #linear interpolation because k=1, ext=0=>extrapolate. Identical to interp1d 

			"""
			#FOR DEBUGGING 
			#visualize interpolation in each bin
			print(f"bin {t}, DM flux = {fluxval[t]} with dloglike={loglike_interp(fluxval[t])}")
			if DMprop[0]==100 and DMprop[1] == 1.3141473626117527e-23:
			#if 1==1:
				print(maxBin)
				xrn =  np.sort(table_norm[t])
				xrn_longer = np.logspace(-25,-9, 100)
				plt.plot(table_norm[t], table_loglike[t], 'o', xrn_longer, loglike_interp(xrn_longer), '-', fluxval[t], loglike_interp(fluxval[t]), 'o')
				plt.xlabel(f'measured flux [ph/cm^2/s/GeV], energy bin {t}')
				plt.ylabel(f'dloglikelihood, energy bin {t}')
				plt.xscale('log')
				plt.show()
				plt.clf()
			"""
	
			LogLike.append(loglike_interp(fluxval[t]))
			del loglike_interp #make sure it clears
		else:
			#print(f"Flux larger than measured: DM flux = {fluxval[t]} compared to measured {max(table_norm[t])}")
			continue

	"""
	#DEBUG
	print(f"dloglike = {np.sum(LogLike)}")
	plt.plot(eref, LogLike, "o")	
	plt.xlabel("Photon energy [GeV]")
	plt.ylabel("dloglike of estimated DM flux")
	plt.title(f"DM Mass = {DMprop[0]} GeV, DM sigmav = {DMprop[1]}")
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	plt.clf()
	"""
			 
	return sum(LogLike)
	
	
############################################################################################################
# MAIN
############################################################################################################
def main(cmd_line):

	'''Here we can set variable for path definitions
		In addition to all desired parameters,
		double check you have the correct LTcube, extdir
	'''

	print("*****running!!!*****")
	
	srcname = cmd_line[1]
	m_BH = float(cmd_line[2]) #black hole mass, units of 10e8 solar masses
	distance = float(cmd_line[3]) #distance to galacy, units of Mpc
	dmMassRange = [int(cmd_line[4]), int(cmd_line[5]), int(cmd_line[6])]
	dmXSRange = [int(cmd_line[7]), int(cmd_line[8]), int(cmd_line[9])]


	sed_suffix = 'free_BG_sed'
	homedir = "/Users/asteinhe/FermiLAT/BHinEGs_DM/run/"

	sedfits =  f"{homedir}{srcname}/output/{srcname}_{sed_suffix}.fits"
	global save_array_path 
	save_array_path = f"{homedir}{srcname}/output/dloglike/{subdir}/"
	makeDir(save_array_path)

	print("Running for {}".format(srcname))

	sed = fits.open(sedfits)[1].data
	refdnde_vec = copy.deepcopy(sed['ref_dnde'])/1.e3 #Differential flux of the reference model evaluated at the bin center (cm−2 s−1 MeV−1 converted to cm-2 s-1 GeV-1)
	#table_normscan1 = copy.deepcopy(sed['norm_scan']) #Array of NxM normalization values for the profile likelihood scan in N energy bins and M scan points. A row-wise multiplication with any of ref columns can be used to convert this matrix to the respective unit.
	table_normscan = copy.deepcopy(sed['norm_scan']) #Array of NxM normalization values for the profile likelihood scan in N energy bins [arrays] and M scan points. A row-wise multiplication with any of ref columns can be used to convert this matrix to the respective unit.
	#convert table_normscan to units of cm^-2 s^-1 GeV^-1
	for t in range(len(refdnde_vec)):
		#table_normscan[t,:] = sed['norm_scan'][t,:]*refdnde_vec[t]
		table_normscan[t] = sed['norm_scan'][t]*refdnde_vec[t]

	table_loglikescan = copy.deepcopy(sed['dloglike_scan']) #Array of NxM delta-loglikelihood values for the profile likelihood scan in N energy bins and M scan points.

	# Convert MeV to GeV
	print("Define eref_vec in GeV...")
	eref_vec = copy.deepcopy(sed['e_ref'])/1.e3 #"reference energy" (= center of bin) of SED energy bins in GeV
	emin_vec = copy.deepcopy(sed['e_min'])/1.e3 #lower edges of SED energy bins in GeV
	emax_vec = copy.deepcopy(sed['e_max'])/1.e3 #upper edges of SED energy bins in GeV
	
	print(eref_vec)
	print(emin_vec)
	
	"""
	#DEBUG
	print(f"ref_dnde = {refdnde_vec}" )
	print(f"e_ref in GeV = {eref_vec}")
	plt.plot(eref_vec, refdnde_vec, "o")
	plt.xlabel("Bin centers [GeV]")
	plt.ylabel("ref_dnde [ph/cm2/s/Gev]")
	plt.yscale('log')
	plt.xscale('log')
	plt.show()
	plt.clf()
	difflux = [refdnde_vec[i]*eref_vec[i]*eref_vec[i] for i in range(len(eref_vec))]
	plt.plot(eref_vec, difflux, "o")
	plt.xlabel("Bin centers [GeV]")
	plt.ylabel("Differential flux (ref_dnde*E*E) [GeV/cm2/s]")
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	plt.clf()
	
	print(f"norm_scan FROM FILE= {table_normscan1}")
	plt.boxplot(table_normscan1.T, showmeans=True)
	plt.xlabel('Energy bin')
	plt.ylabel('Normalization')
	plt.yscale('log')
	plt.show()
	plt.clf()
	
	print(f"table_normscan (converted units)= {table_normscan}")	
	plt.boxplot(table_normscan.T, showmeans=True)
	plt.xlabel('Energy bin')
	plt.ylabel('Normalization*refdnde_vec [ph/cm2/s/GeV]')
	plt.yscale('log')
	plt.show()
	plt.clf()
	
	print(f"dloglike_scan = {table_loglikescan}")
	plt.plot(table_loglikescan.T, "-o")
	plt.xlabel('Scan points')
	plt.ylabel('dloglike per energy bin')
	plt.show()
	plt.clf()
	
	print(f"e_min in GeV = {emin_vec}")
	print(f"e_max in GeV = {emax_vec}")
	print(f"e_ref in GeV = {eref_vec}")
	plt.plot(emin_vec, "o",eref_vec, "o", emax_vec, "o")
	plt.xlabel("Energy bin")
	plt.ylabel("Energy [Gev]")
	plt.yscale('log')
	plt.show()
	plt.clf()
	"""

	#-------------------------------------------
	# Producing Loglike profiling with sigmav.
	#-------------------------------------------

	# create results file
	mass_vec = np.logspace(dmMassRange[0], dmMassRange[1], dmMassRange[2])
	sigmav_vec = np.logspace(dmXSRange[0], dmXSRange[1], dmXSRange[2])
	LogLike_vec = np.zeros(shape=(len(mass_vec),len(sigmav_vec)))

	for u in range(len(mass_vec)):
		for t in range(len(sigmav_vec)):
			#plot energy spectrum with all normalization 
			photSpectrum = diffPhotSpectrum(mass_vec[u], sigmav_vec[t], m_BH, distance, eref_vec) #valid in SED range
			DMprop = [mass_vec[u],sigmav_vec[t], m_BH, distance]
			LogLike_vec[u,t] = funcLogLike(DMprop,photSpectrum,table_normscan,table_loglikescan,emin_vec,emax_vec,eref_vec)

	np.save(save_array_path+'{}_dlike.npy'.format(srcname), LogLike_vec)

################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":

	subdir = "mspPlusDm"
	main(sys.argv)