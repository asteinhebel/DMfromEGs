import os
import numpy as np
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter, SEDPlotter
import matplotlib.pyplot as plt
import matplotlib

#loggerLevel=10 #debug
#loggerLevel=20 #info, default
#loggerLevel=30 #warning

gta = GTAnalysis('config_NGC6251_noSource.yaml')#, loglevel=loggerLevel)
matplotlib.interactive(True)
gta.setup()
gta.write_roi('outputs/fit0_prefit')
gta.print_roi()

#check prefit quality with plots
#define source as point source
resid = gta.residmap('outputs/ngc6251_prefit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
fig = plt.figure(figsize=(14,6))
ROIPlotter(resid['data'],roi=gta.roi).plot(vmin=0,vmax=1000,subplot=121,cmap='magma')
plt.gca().set_title('Data')
ROIPlotter(resid['model'],roi=gta.roi).plot(vmin=0,vmax=1000,subplot=122,cmap='magma')
plt.gca().set_title('Model')
plt.savefig('outputs/plots/prefit_dataModel.png')

#check prefit quality with plots
fig = plt.figure(figsize=(14,6))
ROIPlotter(resid['sigma'],roi=gta.roi).plot(vmin=-5,vmax=5,levels=[-5,-3,3,5],subplot=121,cmap='RdBu_r')
plt.gca().set_title('Significance')
ROIPlotter(resid['excess'],roi=gta.roi).plot(vmin=-100,vmax=100,subplot=122,cmap='RdBu_r')
plt.gca().set_title('Excess')
plt.savefig('outputs/plots/prefit_significanceExcess.png')

#OPTIMIZE ROI MODEL

#This method will iteratively optimize the parameters of all components in the ROI in several stages:
###Simultaneously fitting the normalization of the brightest model components containing at least some fraction of the total model counts (default 95%).
###Individually fitting the normalization of all remaining sources if they have Npred above some threshold (default 1).
###Individually fitting the normalization and shape of any component with TS larger than some threshold (default 25).

#background fit
gta.optimize()
gta.write_roi('outputs/fit1_background')
gta.print_roi()

#check background fit quality with plots
resid = gta.residmap('outputs/ngc6251_postfit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
fig = plt.figure(figsize=(14,6))
ROIPlotter(resid['sigma'],roi=gta.roi).plot(vmin=-5,vmax=5,levels=[-5,-3,3,5],subplot=121,cmap='RdBu_r')
plt.gca().set_title('Significance')
ROIPlotter(resid['excess'],roi=gta.roi).plot(vmin=-100,vmax=100,subplot=122,cmap='RdBu_r')
plt.gca().set_title('Excess')
plt.savefig('outputs/plots/bkgfit_significanceExcess.png')

#check background fit quality with TS map
tsmap_postfit = gta.tsmap('outputs/ngc6251_postfit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
fig = plt.figure(figsize=(14,6))
ROIPlotter(tsmap_postfit['sqrt_ts'],roi=gta.roi).plot(levels=[0,3,5,7],vmin=0,vmax=5,subplot=121,cmap='magma')
plt.gca().set_title('Sqrt(TS)')
ROIPlotter(tsmap_postfit['npred'],roi=gta.roi).plot(vmin=0,vmax=100,subplot=122,cmap='magma')
plt.gca().set_title('NPred')
plt.savefig('outputs/plots/bkgfit_TSNPred.png')


#Find sources outside of the 4FGL catalog
src = gta.find_sources(sqrt_ts_threshold=2)


#Source identified:
#print(gta.roi['PS J2300.3+2913']) #prints TS - don't need to manually identify new sources 
resid = gta.residmap('outputs/ngc6251_newsrcs',model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
fig = plt.figure(figsize=(14,6))
ROIPlotter(resid['sigma'],roi=gta.roi).plot(vmin=-5,vmax=5,levels=[-5,-3,3,5],subplot=121,cmap='RdBu_r')
plt.gca().set_title('Significance')
ROIPlotter(resid['excess'],roi=gta.roi).plot(vmin=-100,vmax=100,subplot=122,cmap='RdBu_r')
plt.gca().set_title('Excess')
plt.savefig('outputs/plots/newsrcs_significanceExcess.png')

tsmap_newsrcs = gta.tsmap('outputs/ngc6251_newsrcs',model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
fig = plt.figure(figsize=(14,6))
ROIPlotter(tsmap_newsrcs['sqrt_ts'],roi=gta.roi).plot(levels=[0,3,5,7],vmin=0,vmax=5,subplot=121,cmap='magma')
plt.gca().set_title('Sqrt(TS)')
ROIPlotter(tsmap_newsrcs['npred'],roi=gta.roi).plot(vmin=0,vmax=100,subplot=122,cmap='magma')
plt.gca().set_title('NPred')
plt.savefig('outputs/plots/newsrcs_TSNPred.png')

#SPECTRAL ANALYSIS OF SOURCE

#Free parameters of source and everything within 3deg
gta.free_sources(distance=3.0,pars='norm')
gta.free_sources(distance=3.0,pars='shape',minmax_ts=[100.,None])
#refit
fit_results = gta.fit()
gta.write_roi('outputs/fit3_sourceFree')
gta.print_roi()
print(f"dloglike: {fit_results['dloglike']}")


#execute spectral analysis - bin-by-bin likelihood scan
sed_ngc6251 = gta.sed('ngc6251')
sed_j1540 = gta.sed('4FGL J1540.1+8155') #for comparison
gta.write_roi('outputs/fit_sed')

#inspect fit results
# E^2 x Differential flux ULs in each bin in units of MeV cm^{-2} s^{-1}
print(sed_ngc6251['e2dnde_ul95'])
#dloglike_scan array saves likelihood in each bin
e2dnde_scan = sed_ngc6251['norm_scan']*sed_ngc6251['ref_e2dnde'][:,None]

for b in range(len(sed_ngc6251['dloglike_scan'])):
	plt.figure()
	plt.plot(e2dnde_scan[b],sed_ngc6251['dloglike_scan'][b]-np.max(sed_ngc6251['dloglike_scan'][b]),label='Log likelihood')
	plt.gca().set_ylim(-5,1)
	plt.gca().axvline(sed_ngc6251['e2dnde_ul95'][b],color='k',label='Flux UL')#flux UL  in this bin
	plt.gca().axhline(-2.71/2.,color='r',label='Reduction of -2.71/2 from max')
	plt.xlabel('Energy [MeV]')
	plt.ylabel('Delta Log likelihood')
	plt.title('Bin '+str(b)+' of '+str(len(sed_ngc6251['dloglike_scan'])))
	plt.legend(loc='best')
	plt.savefig('outputs/plots/likelihood_bin'+str(b)+'.png')
	plt.close()

#visualize whole SED
fig = plt.figure(figsize=(14,4))
ylim=[1E-8,1E-5]
xlim=[1E3,1E6]
fig.add_subplot(121)
SEDPlotter(sed_ngc6251).plot()
plt.gca().set_ylim(ylim) 

fig.add_subplot(122)
SEDPlotter(sed_j1540).plot()
plt.gca().set_ylim(ylim) 
plt.savefig('outputs/plots/likelihood_fullSED_ylim.png')


fig = plt.figure(figsize=(14,4))

fig.add_subplot(121)
SEDPlotter(sed_ngc6251).plot(showlnl=True,ylim=ylim)
plt.gca().set_ylim(ylim)

fig.add_subplot(122)
SEDPlotter(sed_j1540).plot(showlnl=True,ylim=ylim)
plt.gca().set_ylim(ylim)
plt.savefig('outputs/plots/likelihood_fullSED_deltaLogLike.png')


#SET UPPER LIMITS

import pyLikelihood

# Load the sed data structure
data = sed_ngc6251

# Instantiate a DM Fit Function for a DM particle spectrum given the following parameters
### must define model
# Mass = 100 GeV
# Cross-Section: 3 x 10^{-26} cm^{3} s^{-1}
# J-Factor: 10^19 GeV^2 cm^{-5}
# Channel: b-bbar
dmf = pyLikelihood.DMFitFunction()
dmf.readFunction(os.path.expandvars('gammamc_dif.dat'))
dmf.setParam('norm',1E19)
dmf.setParam('sigmav',3E-26)
dmf.setParam('mass',100.0)
dmf.setParam('bratio',1.0)
dmf.setParam('channel0',4)

def integrate_eflux(fn,ebins,nstep=10):
    #Compute energy flux within a sequence of energy bins.
    
    loge = np.linspace(ebins[0],ebins[-1],100)
    dfde = [fn(pyLikelihood.dArg(10**x)) for x in loge]        
    dfde = np.array(dfde)
    x = ebins
    dx = (x[1:] - x[:-1])

    yedge = x[1:,np.newaxis] + np.linspace(0,1,nstep)[np.newaxis,:]*dx[:,np.newaxis] 
    dy = 10**yedge[:,1:]-10**yedge[:,:-1]
    y = 0.5*(yedge[:,1:]+yedge[:,:-1])
    eflux = np.interp(np.ravel(y),loge,dfde)
    eflux = np.sum(eflux.reshape(y.shape)*10**y*dy,axis=1)

    return eflux

class SEDLike(object):

    def __init__(self,sed):
        self._sed = sed
        self._eflux_scan = sed['norm_scan']*sed['ref_eflux'][:,None]

    def __call__(self,eflux):
        lnl = np.zeros(eflux.shape)
        for i, ectr in enumerate(self._sed['e_ctr']):
            v = np.interp(eflux[i],
                          self._eflux_scan[i],
                          self._sed['dloglike_scan'][i])
            lnl[i] += v
        return np.sum(lnl,axis=0)

ebins = np.log10(np.array(list(data['e_min']) + list([data['e_max'][-1]])))
eflux = integrate_eflux(dmf,ebins)
sigmav = 3.E-26*np.logspace(-3.,1.,101)
eflux = eflux[:,np.newaxis]*np.logspace(-3,1,101)[np.newaxis,:]

slike = SEDLike(data)
lnl = slike(eflux)
lnl -= np.max(lnl)

# Plot log-likelihood profile

plt.figure()
plt.plot(sigmav,lnl)
plt.gca().set_xscale('log')
plt.gca().axhline(-2.71/2.,color='k')
plt.gca().set_ylim(-4,1)
plt.gca().set_xlabel('Cross Section [cm$^{3}$ s$^{-1}$]')

sigmav_ul = float(np.interp(2.71/2.,-lnl,sigmav))

print(f'Sigma-V Upper Limit: {sigmav_ul}')
plt.savefig('outputs/plots/dmLimits_defaultModel.png')
