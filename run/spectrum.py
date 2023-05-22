import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def convertGeVErg(g):
	erg = 0.0000000000016 * g * 1000000. #g in GeV but conversion for just eV
	return erg

def diffEn(m,e):
	E = convertGeVErg(e)
	frac = (E/m)**2
	return frac * (E/m) * np.exp(-frac)
	
############################################################################################################
#MAIN
############################################################################################################
def main():

	#mass_vec = [10]
	mass_vec =  [1,5,10,25,50,100,250,500,1000,5000,10000] #GeV
	gam_vec =  np.logspace(-1,4,50) #GeV 
	interp_x =  np.logspace(-1,4,200) #GeV 

	specMatrix = np.zeros((len(mass_vec), len(gam_vec)))
	interp = []

	for i,m in enumerate(mass_vec):
		for j,e in enumerate(gam_vec):
			specMatrix[i][j] = diffEn(m,e)
		interp.append(np.interp(interp_x, gam_vec, specMatrix[i]))
		plt.plot(interp_x, interp[-1], label=f"{m:.0f} GeV")
	
	plt.legend(loc="best")
	plt.xscale('log')
	plt.yscale('log')
	plt.title("Differential Energy Flux")
	plt.xlabel("photon energy [GeV]")
	plt.ylabel("estimated DM flux [erg/cm^2/s/GeV]")
	plt.show()
	plt.clf()
	
	#make finer granularity
	photMatrix = np.zeros((len(mass_vec), len(interp_x)))
	for i,m in enumerate(mass_vec):
		for j,e in enumerate(interp_x):
			photMatrix[i][j] = interp[i][j]/e
		plt.plot(interp_x, photMatrix[i], '-', label=f"{m:.2f} GeV")
	
	plt.legend(loc="best")
	plt.xscale('log')
	plt.yscale('log')
	plt.title("Differential Photon Flux")
	plt.xlabel("photon energy [GeV]")
	plt.ylabel("estimated DM flux [/cm^2/s]")
	plt.show()

################################################################################################
################################################################################################
if __name__=="__main__":
	
	main()