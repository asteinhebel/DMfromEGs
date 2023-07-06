# -*- coding: utf-8 -*-
#imports
from time import sleep
from random import randint
import resource
import random,shutil,yaml
import os,sys
from math import *
import numpy as np
import pandas as pd
import matplotlib as matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter, SEDPlotter

print("create_vectors")
obsEnergies = np.logspace(-2,3,500) #GeV
mass_vec =  np.logspace(-1,3,5) #GeV
sigmav_vec = np.logspace(-28,-22,5) #cm^3/s


def diffEn(m,e):
    frac = (e/m)**2
    flux = frac * (e/m) * np.exp(-frac)
    return flux

def full_diffPhotSpectrum(m, sv, mBH, dist):
    #Calculate dE/dN of DM spectrum, units of ph/cm^2/s/GeV (after multiplication with DM sigma), return as array of points that span LAT sensitivity
    print(mBH, sv, dist, m)
    norm = 2.e-20 * pow(mBH, (5./2.)) * (sv/1e-26) * pow(1./m, 2.) * pow(1./dist,2.) * (2 / m) # ph/cm2/s/GeV
    spectrum = np.array([norm * diffEn(m, en) for en in obsEnergies])
    
    return spectrum
    
def generate_flux_DM_model_file(DMprop, spec, folder):

    dict = {'Energy[MeV]'           : obsEnergies*1000,
            'DiffFlux[ph/cm2/s/GeV]': spec
            }
            
    df = pd.DataFrame.from_dict(dict,)
    df.to_csv('%s/Spec_DM_m%.1f_sv%.1f.txt'%(folder, DMprop[0], -np.log10(DMprop[1])),  header=False, index=False, sep=' ')
    return 0


def create_DMspectra_files(srcname, mBH, distance):
    '''
    '''
    print("Creating the files")
    
    sed_suffix = 'free_BG_sed'
    homedir = "/Users/mnegro/MyDocuments/_FERMI/DM_EllipticalGalaxies/DMfromEGs/run/"
    
    src_folder = f"{homedir}{srcname}/"
    sedfits =  f"{homedir}{srcname}/output/{srcname}_{sed_suffix}.fits"
    
    for u in range(len(mass_vec)):
    
        for t in range(len(sigmav_vec)):
            print('DM mass, DM xsection: ', mass_vec[u], sigmav_vec[t])
            photSpectrum = full_diffPhotSpectrum(mass_vec[u], sigmav_vec[t], mBH, distance)
            DMprop = [mass_vec[u], sigmav_vec[t], mBH, distance]
            generate_flux_DM_model_file(DMprop, photSpectrum, src_folder)

def ang_sep(ra0, dec0, ra1, dec1):
    C = np.pi/180.
    d0 = C * dec0
    d1 = C * dec1
    r12 = C * (ra0 - ra1)
    cd0 = np.cos(d0)
    sd0 = np.sin(d0)
    cd1 = np.cos(d1)
    sd1 = np.sin(d1)
    cr12 = np.cos(r12)
    sr12 = np.sin(r12)
    num = np.sqrt((cd0 * sr12) ** 2 + (cd1 * sd0 - sd1 * cd0 * cr12) ** 2)
    den = sd0 * sd1 + cd0 * cd1 * cr12
    return np.arctan2(num, den) / C

def savePlots(gta,name,svdir):
	resid = gta.residmap(name,model={'SpatialModel' : 'PointSource', 'Index' : 2.0})
	fig = plt.figure(figsize=(14,6))
	ROIPlotter(resid['data'],roi=gta.roi).plot(vmin=0,vmax=1000,subplot=121,cmap='magma')
	plt.gca().set_title('Data')
	ROIPlotter(resid['model'],roi=gta.roi).plot(vmin=0,vmax=1000,subplot=122,cmap='magma')
	plt.gca().set_title('Model')
	plt.savefig(svdir+name+'_map.png')

	#check prefit quality with plots
	fig = plt.figure(figsize=(14,6))
	ROIPlotter(resid['sigma'],roi=gta.roi).plot(vmin=-5,vmax=5,levels=[-5,-3,3,5],subplot=121,cmap='RdBu_r')
	plt.gca().set_title('Significance')
	ROIPlotter(resid['excess'],roi=gta.roi).plot(vmin=-100,vmax=100,subplot=122,cmap='RdBu_r')
	plt.gca().set_title('Excess')
	plt.savefig(svdir+name+'_significanceExcess.png')

def makeDir(dirpath, overwrite=False):
    if(os.path.isdir(dirpath)==True):
        print(f"Path to {dirpath} exists")
        if overwrite:
            shutil.rmtree(dirpath)
            os.system('mkdir -p %s' %dirpath)
        else:
            pass
    else:
        print(f"path to {dirpath} does not exist, creating...")
        os.system('mkdir -p %s' %dirpath)

############################################################################################################
#MAIN
############################################################################################################
def main(cmd_line):

    '''Here we can set variable for path definitions
       In addition to all desired parameters, 
       double check you have the correct LTcube, extdir
    '''
    
    homedir = '/Users/mnegro/MyDocuments/_FERMI/DM_EllipticalGalaxies/DMfromEGs/'

####################################
	#Get inputs
####################################
    print("*****running!!!*****")
    print(cmd_line)
    srcname = cmd_line[1]
    ra      = cmd_line[2]
    dec     = cmd_line[3]
    m_BH    = cmd_line[4]
    dist_bh = cmd_line[5]

    outdir = homedir + 'run/' + srcname + '/'
    makeDir(outdir, overwrite=False)
    os.chdir(outdir)

    #make sure all outputs exist
    outputdir = homedir + 'run/' +srcname + '/output/'
    makeDir(outputdir, overwrite=False)
    plotDir = homedir + 'run/' +srcname + '/output/plots/'
    makeDir(plotDir, overwrite=False)

#
#####################################
#    #Create config file
#####################################
#    create_DMspectra_files(srcname, m_BH, dist_bh)
#    print('---> done creating spectrum files...\n\n\n')
#
#    ft1    = outdir + 'ft1.txt' #data files
#    ft2    = homedir + 'data/ft2_sc.fits' #SC file
#    ltcube = homedir + 'data/ltcube_00.fits'
#    galdiff = '/Users/mnegro/opt/anaconda3/envs/fermi/share/fermitools/refdata/fermi/galdiffuse/gll_iem_v07.fits'
#    isodiff = '/Users/mnegro/opt/anaconda3/envs/fermi/share/fermitools/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt'
#
#    with open('%s.yaml' % srcname, 'w') as yml: #create the yml
#        yml.write("logging:\n")
#        yml.write("  verbosity : 3\n")
#        yml.write("  chatter : 3\n")
#        yml.write("#--------#\n")
#        yml.write("fileio:\n")
#        yml.write("  outdir : output\n")
#        yml.write("  logfile : %s\n" % srcname)
#        yml.write("  usescratch : False\n")
#        yml.write("  scratchdir : scratch\n")
#        yml.write("#--------#\n")
#        yml.write("data:\n")
#        yml.write("  evfile : '%s'\n" % ft1)
#        yml.write("  scfile : '%s'\n" % ft2)
#        yml.write("  ltcube : '%s'\n" % ltcube)
#        yml.write("#--------#\n")
#        yml.write("binning:\n")
#        yml.write("  roiwidth : 10\n")
#        yml.write("  binsz : 0.08\n")
#        yml.write("  binsperdec : 8\n")
#        yml.write("#--------#\n")
#        yml.write("selection:\n")
#        yml.write("  emin : 1000\n")
#        yml.write("  emax : 1000000\n")
#        yml.write("  zmax : 100\n")
#        yml.write("  ra   : %.3f\n"%ra)
#        yml.write("  dec  : %.3f\n"%dec)
##        yml.write("  target : '%s'\n" % srcname)
#        yml.write("  radius : 15\n")
#        yml.write("  tmin : 239557417\n")
#        yml.write("  tmax : 665410983\n")
#        yml.write("  evclass : 128\n")
#        yml.write("  evtype : 3\n")
#        yml.write("  filter : 'DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52'\n")
#        yml.write("#--------#\n")
#        yml.write("gtlike:\n")
#        yml.write("  edisp : True\n")
#        yml.write("  edisp_disable : ['isodiff']\n")
#        yml.write("  irfs : 'P8R3_SOURCE_V3'\n")
#        yml.write("#--------#\n")
#        yml.write("model:\n")
#        yml.write("  src_radius : 15\n")
#        yml.write("  src_roiwidth : 15\n")
#        yml.write("  galdiff : '%s'\n" %galdiff)
#        yml.write("  isodiff : '%s'\n" %isodiff)
#        yml.write("  catalogs :\n")
#        yml.write("    - '4FGL-DR3'\n")
#        yml.write("#--------#\n")
#        yml.write("components: null\n")
#        yml.write("plotting:\n")
#        yml.write("  format : png\n")
#        yml.write("#--------#\n")
#        yml.write("sed:\n")
#        yml.write("  use_local_index : True")
#    yml.close()
##
#####################################
#    #Configure/run GTAnalysis
#####################################
    gta = GTAnalysis('%s.yaml' % srcname, logging={'verbosity' : 3})
    gta.setup()
#

#    gta.optimize()
#    gta.print_roi()
#    gta.free_source('galdiff')
#    gta.free_source('isodiff')
#    gta.free_sources(minmax_ts=[25,None],pars='norm', distance=5.0)
#    gta.free_sources(minmax_ts=[500,None],distance=7.0)
#    gta.print_roi()
#    gta.fit()
#    gta.print_roi()
#    roi_base = gta.write_roi('fit_model_base')
#    savePlots(gta, 'fit_model_base', plotDir)
#
#    print("***** loading the base model !!!   *****")
#    gta.load_roi('fit_model_base')
#    gta.print_roi()
    #This finds new sources with TS>25
#    gta.free_sources(free=False)
#    sdict = gta.find_sources(sqrt_ts_threshold=2.0,min_separation=0.15,tsmap_fitter='tsmap')
#    gta.free_sources(free=True)
#    gta.free_sources(minmax_ts=[25, None], free=False, pars='norm')
#    gta.print_model()
#    gta.fit(optimizer='NEWMINUIT')
#    gta.free_sources(free=False)
#    gta.print_model()
#    savePlots(gta,'fit_model_base_addsrcTS25', plotDir)
#    gta.write_roi('fit_model_base_addsrcTS25', make_plots=True, save_model_map=True)
#    os.system('mv %s*.png %s.'%( homedir+'run/', plotDir))
    
    print("***** loading the base model + new src (TS > 4) !!!   *****")
    gta.load_roi('fit_model_base_addsrcTS25')
    gta.print_roi()
    
    _dll_ = []
    for i, m in enumerate(mass_vec):
        dll_ = []
        for j, sv in enumerate(sigmav_vec):
            print("***** loading the base model + new srcs!!!   *****")
            gta.load_roi('fit_model_base_addsrcTS25')
            gta.free_sources(free=False)

            spec_filename = homedir + 'run/' + srcname + '/Spec_DM_m%.1f_sv%.1f.txt'%(m, -np.log10(sv))
            src_dict = { 'ra' : ra, 'dec' : dec, 'eflux': 1e-14,
                         'SpectrumType' : 'FileFunction',
                         'Spectrum_Filename': spec_filename
                        }
            gta.add_source(srcname + '_%i%i'%(i, j), src_dict, free=True, init_source=True, save_source_maps=True)
#            gta.optimize()
#            gta.free_sources(minmax_ts=[None, 9], free=True)
#            gta.free_source('galdiff')
            gta.free_source('isodiff')
            gta.print_model()
            
            fit_results = gta.fit(optimizer='NEWMINUIT')
            gta.print_model()
            print('-DeltaLogLike', fit_results['dloglike'])
            dll_.append(fit_results['dloglike'])

            gta.delete_source(srcname+'_%i%i'%(i, j))
            gta.print_model()
            
        _dll_.append(dll_)
    _dll_ = np.array(_dll_)
    _TS_ = -2 * _dll_
    
    plt.imshow(_TS_, origin='lower', aspect='auto' )
    plt.savefig('AAAAA.png', dpi=200)
#            savePlots(gta,'fit_model_base_addsrcTS25', plotDir)
#            gta.write_roi('fit_model_base_addsrcTS25', make_plots=True, save_model_map=True)
    plt.show()

#    ##########
#    #Model 1
#    gta.optimize()
#    gta.print_roi()
#    gta.free_source('galdiff')
#    gta.free_source('isodiff')
#    gta.free_sources(minmax_ts=[25,None],pars='norm', distance=5.0)
#    gta.free_sources(minmax_ts=[500,None],distance=7.0)
#    gta.free_source(srcname)
#    gta.fit()
#    gta.print_roi()
#    gta.write_roi('fit_model_1')
#    savePlots(gta,'fit_model_1', plotDir)
#
#    ##########
#    #Model 2
#    gta.delete_source(srcname)
#    model2 = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
#    finder2 = gta.find_sources(prefix='find2',model=model2,sqrt_ts_threshold=4.0,
#                min_separation=0.5,max_iter=10,sources_per_iter=20)
#                #,tsmap_fitter='tsmap')
#    gta.print_roi()
#    gta.write_roi('fit_model_2')
#    #savePlots(gta,'fit_model_2', plotDir)
#
#    ##########
#    #Model 3
#    names2    = np.array([finder2['sources'][i]['name'] for i in range(len(finder2['sources']))])
#    ra_ns2    = np.array([finder2['sources'][i]['ra'] for i in range(len(finder2['sources']))])
#    dec_ns2   = np.array([finder2['sources'][i]['dec'] for i in range(len(finder2['sources']))])
#    r95_ns2   = np.array([finder2['sources'][i]['pos_r95'] for i in range(len(finder2['sources']))])
#    dist_sep = '%.2f' %0
#    namee=''
#
#    #check ang separation if more than one source is found close to ra and dec
#    if len(names2)>0:
#        for i in range(len(names2)):
#            sepp2=ang_sep(ra_ns2[i],dec_ns2[i],ra,dec)
#            if sepp2 < r95_ns2[i] and sepp2 < 0.2:
#                print(names2[i],sepp2,r95_ns2[i])
#                dist_sep = '%.2f' %sepp2
#                namee = names2[i]
#    if namee:
#        srcname=namee
#    else:
#        gta.add_source(srcname,{ 'ra' : ra, 'dec' : dec,
#            'SpectrumType' : 'PowerLaw', 'Index' : 2.0,
#            'Scale' : 1000, 'Prefactor' : 1e-11,
#            'SpatialModel' : 'PointSource' })
#
#    gta.optimize()
#    gta.print_roi()
#    gta.free_source('galdiff')
#    gta.free_sources(minmax_ts=[500,None],distance=7.0)
#    gta.free_source(srcname)
#    gta.fit()
#    gta.print_roi()
#    gta.write_roi('fit_model_3')
#    savePlots(gta,'fit_model_3',plotDir)
#
#####################################
#    #Get params from model 3 and save output
#####################################
#    p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
#    src = p['sources'][srcname]
#    if str(src['ts'])=='nan':   # To check non-convergence
#        print('****************************')
#        print('Fit has not converged')
#        print('****************************')
#        gta.free_sources(minmax_ts=[None,100],free=False)
#        gta.free_source(srcname)
#        gta.fit(tol=1e-8)
#        gta.write_roi('fit_model_3')
#        p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
#        rc = p['sources'][srcname]
#    else:
#        print('****************************')
#        print('Fit has converged')
#        print('****************************')
#
#    """
#    if src['ts'] > 15 and np.fabs(src['param_values'][1]) > 3.0:
#        model4 = {'Index' : 2.8, 'SpatialModel' : 'PointSource'}
#        print('****************************')
#        print('Source still bad TS={} Index={}?'.format(src['ts'],src['param_values'][1]))
#        print('****************************')
#        gta.delete_source(srcname)
#        mapp=gta.tsmap('fit_no_source_final',model=model4)
#        finder4 = gta.find_sources(prefix='find4',model=model4,sqrt_ts_threshold=4.0,
#                      min_separation=1.0,max_iter=10,sources_per_iter=20,
#                      tsmap_fitter='tsmap')
#        gta.add_source(srcname,{ 'ra' : ra, 'dec' : dec,
#                    'SpectrumType' : 'PowerLaw', 'Index' : 2.0,
#                    'Scale' : 1000, 'Prefactor' : 1e-11,
#                    'SpatialModel' : 'PointSource' })
#        gta.free_sources(minmax_ts=[100,None],pars='norm',distance=5.0)
#        gta.free_sources(minmax_ts=[200,None],distance=7.0)
#        gta.free_source(srcname)
#        gta.fit()
#        gta.write_roi('fit_model_3')
#        p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
#        src = p['sources'][srcname]
#    """
#
#    TS = '%.2f' %src['ts']
#    Flux = '%.2e' %src['eflux']
#    Flux_err = '%.2e' %src['eflux_err']
#    Flux_UL = '%.2e' %src['eflux_ul95']
#    Index = '%.2f' %(np.fabs(src['param_values'][1]))
#    Index_err = '%.2f' %src['param_errors'][1]
#    f = open('%s_Param.txt' % srcname, 'w')
#    f.write(str(srcname) + "\t" + str(dist_sep) + "\t" + str(Flux) + "\t" + str(Flux_err) + "\t" + str(Index) + "\t" + str(Index_err) + "\t" + str(Flux_UL) + "\t" + str(TS) + "\n")
#    f.close()
    
    return
    
################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":

    targets = ["NGC4889", "NGC4649", "NGC1407"]
    ra      = [195.034, 190.917, 55.0496    ]
    dec     = [27.977, 11.5526, -18.5804    ]
    m_bh    = [208.0, 47.3, 46.5] #1e8 solar masses
    distance= [102.0, 16.5, 29.0] #Mpc
    

    main(['placeholder', targets[0], ra[0], dec[0], m_bh[0], distance[0]])
