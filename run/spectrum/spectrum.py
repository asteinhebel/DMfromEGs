import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#Define test galaxy - NGC4649
mBH = 47.3
dist = 16.5
sigmav = 1e-26

def convertGeVErg(g):
	erg = 0.000000000001602 * g * 1000000000. #g in GeV but conversion for just eV
	return erg
	
def normalization_flux(DMmass):    
    sig26 = sigmav * 1e26 #scale to get sigma ^-26
    norm = (2.e-20)*pow(mBH,5./2.)*sig26*pow(1./DMmass,2.)*pow(1./dist,2.)	
    return norm

def diffEn(m,e):
	frac = (e/m)**2
	norm = 2 * normalization_flux(m) / m 
	flux = norm * frac * (e/m) * np.exp(-frac)
	return e/m, flux
	
############################################################################################################
#MAIN
############################################################################################################
def main():

	#Define placeholder vectors
	mass_vec =  [1,5,10,25,50,100,250,500,1000,5000,10000] #GeV
	gam_vec =  np.logspace(-1,5,50) #GeV 
	interp_x =  np.logspace(-1,5,200) #GeV 

	specMatrix = np.zeros((len(mass_vec), len(gam_vec)))
	emMatrix = np.zeros((len(mass_vec), len(gam_vec)))
	interp = []
	
	#Derive/plot differential spectrum
	for i,m in enumerate(mass_vec):
		for j,e in enumerate(gam_vec):
			emMatrix[i][j],specMatrix[i][j] = diffEn(m,e)
		interp.append(np.interp(interp_x, gam_vec, specMatrix[i]))
		plt.plot(interp_x, interp[-1], label=f"{m:.0f} GeV")
	
	plt.legend(loc="best")
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim([10e-40, 10e-18])
	plt.title("Differential Photon Flux, dN/dE")
	plt.xlabel("photon energy [GeV]")
	plt.ylabel("estimated DM flux [ph/cm^2/s/GeV]")
	#plt.show()
	plt.savefig("ngc4649_spectrum.png")
	plt.clf()
	
	#plot ratio plot like in jeremy's paper : E*F_E (=dN/dE) vs E/m 
	for p in range(len(emMatrix)):
		plt.plot(emMatrix[p], specMatrix[p], label=f"{mass_vec[p]} GeV")
		
	plt.legend(loc="best")
	plt.xscale('log')
	plt.xlim([0.1,10])
	plt.yscale('log')
	plt.ylim([10e-30,10e-18])
	plt.title("Scaled Differential Flux")
	plt.xlabel("E/m")
	plt.ylabel("estimated DM flux [ph/cm^2/s/GeV]")
	#plt.show()
	plt.savefig("ngc4649_spectrum_Em.png")
	plt.clf()
	

################################################################################################
################################################################################################
if __name__=="__main__":
	
	main()