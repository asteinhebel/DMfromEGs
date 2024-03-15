# -*- coding: utf-8 -*-
#imports
from time import sleep
from random import randint
import resource
import random,shutil,yaml
import os,sys
from math import *
import numpy as np
import matplotlib as matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter, SEDPlotter


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

def makeDir(dirpath):
    if(os.path.isdir(dirpath)==True):
        print(f"Path to {dirpath} exists")
        shutil.rmtree(dirpath)
        os.system('mkdir -p %s' %dirpath) 
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
    
    homedir = '/Users/asteinhe/FermiLAT/BHinEGs_DM/'

####################################
	#Get inputs
####################################
    print("*****running!!!*****")
    srcname = cmd_line[1]
    ra      = float(cmd_line[2])
    dec     = float(cmd_line[3])
    
    outdir=homedir+'run/'+srcname+'/'
    makeDir(outdir)
    os.chdir(outdir)
    
    #make sure all outputs exist   
    makeDir(homedir+'run/'+srcname+'/output/')
    makeDir(homedir+'run/'+srcname+'/output/plots/')
    plotDir = homedir+'run/'+srcname+'/output/plots/'

####################################
    #Create config file
####################################
    ft1 = homedir+'data/'+srcname+'/ft1.fits' #data files
    ft2 = homedir+'data/ft2_sc.fits' #SC file
    ltcube = homedir+'data/ltcube_00.fits' 

    galdiff = '/Users/asteinhe/opt/anaconda3/envs/fermi_test/share/fermitools/refdata/fermi/galdiffuse/gll_iem_v07.fits'
    isodiff = '/Users/asteinhe/opt/anaconda3/envs/fermi_test/share/fermitools/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt'  
    
    with open('%s.yaml' % srcname, 'w') as yml: #create the yml
        yml.write("logging:\n")
        yml.write("  verbosity : 3\n")
        yml.write("  chatter : 3\n")
        yml.write("#--------#\n")
        yml.write("fileio:\n")
        yml.write("  outdir : output\n")
        yml.write("  logfile : %s\n" % srcname)
        yml.write("  usescratch : False\n")
        yml.write("  scratchdir : scratch\n")
        yml.write("#--------#\n")
        yml.write("data:\n")
        yml.write("  evfile : '%s'\n" % ft1)
        yml.write("  scfile : '%s'\n" % ft2)
        yml.write("  ltcube : '%s'\n" % ltcube)
        yml.write("#--------#\n")
        yml.write("binning:\n")
        yml.write("  roiwidth : 10\n")
        yml.write("  binsz : 0.08\n")
        yml.write("  binsperdec : 4\n") #8
        yml.write("#--------#\n")
        yml.write("selection:\n")
        yml.write("  emin : 500\n") 
        yml.write("  emax : 1000000\n")
        yml.write("  zmax : 100\n")
        yml.write("  target : '%s'\n" % srcname)
        yml.write("  radius : 15\n")
        yml.write("  tmin : 239557417\n")
        yml.write("  tmax : 712597417\n") 
        yml.write("  evclass : 128\n")
        yml.write("  evtype : 3\n")
        yml.write("  filter : 'DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52'\n")
        yml.write("#--------#\n")
        yml.write("gtlike:\n")
        yml.write("  edisp : True\n")
        yml.write("  edisp_disable : ['isodiff']\n")
        yml.write("  irfs : 'P8R3_SOURCE_V3'\n")
        yml.write("#--------#\n")
        yml.write("model:\n")
        yml.write("  src_radius : 15\n")
        yml.write("  src_roiwidth : 15\n")
        yml.write("  galdiff : '%s'\n" %galdiff)
        yml.write("  isodiff : '%s'\n" %isodiff)
        yml.write("  catalogs :\n")
        yml.write("    - '4FGL-DR3'\n")
        yml.write("  sources :\n")
        yml.write("    - { 'name' : '%s', 'ra' : %s, 'dec' : %s, 'SpectrumType' : PowerLaw }\n" % (srcname,ra,dec)) #'Spectrum_Filename', '<txtFile>, 'SpatialModel','PointSource'
        yml.write("#--------#\n")
        yml.write("components: null\n")
        yml.write("plotting:\n")
        yml.write("  format : png\n")
        yml.write("#--------#\n")
        yml.write("sed:\n")
        yml.write("  use_local_index : True")
    yml.close()
    
####################################
    #Configure/run GTAnalysis
####################################
    gta = GTAnalysis('%s.yaml' % srcname,logging={'verbosity' : 3})
    gta.setup()
    
    ##########
    #Model 1
    gta.optimize()
    gta.print_roi()
    gta.free_source('galdiff')
    gta.free_source('isodiff')
    gta.free_sources(minmax_ts=[25,None],pars='norm', distance=5.0)
    gta.free_sources(minmax_ts=[500,None],distance=7.0)
    gta.free_source(srcname)
    gta.fit()
    gta.print_roi()
    gta.write_roi('fit_model_1')
    savePlots(gta,'fit_model_1', plotDir)
    
    ##########
    #Model 2
    gta.delete_source(srcname)
    model2 = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    finder2 = gta.find_sources(prefix='find2',model=model2,sqrt_ts_threshold=4.0,
                min_separation=0.5,max_iter=10,sources_per_iter=20)
                #,tsmap_fitter='tsmap')
    gta.print_roi()
    gta.write_roi('fit_model_2')
    #savePlots(gta,'fit_model_2', plotDir)

    ##########
    #Model 3
    names2    = np.array([finder2['sources'][i]['name'] for i in range(len(finder2['sources']))])
    ra_ns2    = np.array([finder2['sources'][i]['ra'] for i in range(len(finder2['sources']))])
    dec_ns2   = np.array([finder2['sources'][i]['dec'] for i in range(len(finder2['sources']))])
    r95_ns2   = np.array([finder2['sources'][i]['pos_r95'] for i in range(len(finder2['sources']))])
    dist_sep = '%.2f' %0
    namee=''
    
    #check ang separation if more than one source is found close to ra and dec
    if len(names2)>0:
        for i in range(len(names2)):
            sepp2=ang_sep(ra_ns2[i],dec_ns2[i],ra,dec)
            if sepp2 < r95_ns2[i] and sepp2 < 0.2:
                print(names2[i],sepp2,r95_ns2[i])
                dist_sep = '%.2f' %sepp2
                namee = names2[i]
    if namee:
        srcname=namee
    else:
        gta.add_source(srcname,{ 'ra' : ra, 'dec' : dec,
            'SpectrumType' : 'PowerLaw', 'Index' : 2.0,
            'Scale' : 1000, 'Prefactor' : 1e-11,
            'SpatialModel' : 'PointSource' })
            
    gta.optimize()
    gta.print_roi()
    gta.free_source('galdiff')
    gta.free_sources(minmax_ts=[500,None],distance=7.0)
    gta.free_source(srcname)
    gta.fit()
    gta.print_roi()
    gta.write_roi('fit_model_3')
    savePlots(gta,'fit_model_3',plotDir)
    
####################################
    #Get params from model 3 and save output
####################################  
    p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
    src = p['sources'][srcname]
    if str(src['ts'])=='nan':   # To check non-convergence
        print('****************************')
        print('Fit has not converged')
        print('****************************')
        gta.free_sources(minmax_ts=[None,100],free=False)
        gta.free_source(srcname)
        gta.fit(tol=1e-8)
        gta.write_roi('fit_model_3')
        p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
        rc = p['sources'][srcname]
    else:
        print('****************************')
        print('Fit has converged')
        print('****************************')

    """
    if src['ts'] > 15 and np.fabs(src['param_values'][1]) > 3.0:
        model4 = {'Index' : 2.8, 'SpatialModel' : 'PointSource'}
        print('****************************')
        print('Source still bad TS={} Index={}?'.format(src['ts'],src['param_values'][1]))
        print('****************************')
        gta.delete_source(srcname)
        mapp=gta.tsmap('fit_no_source_final',model=model4)
        finder4 = gta.find_sources(prefix='find4',model=model4,sqrt_ts_threshold=4.0,
                      min_separation=1.0,max_iter=10,sources_per_iter=20,
                      tsmap_fitter='tsmap')
        gta.add_source(srcname,{ 'ra' : ra, 'dec' : dec,
                    'SpectrumType' : 'PowerLaw', 'Index' : 2.0,
                    'Scale' : 1000, 'Prefactor' : 1e-11,
                    'SpatialModel' : 'PointSource' })
        gta.free_sources(minmax_ts=[100,None],pars='norm',distance=5.0)
        gta.free_sources(minmax_ts=[200,None],distance=7.0)
        gta.free_source(srcname)
        gta.fit()
        gta.write_roi('fit_model_3')
        p = np.load('output/fit_model_3.npy', allow_pickle=True).flat[0]
        src = p['sources'][srcname]
    """

    TS = '%.2f' %src['ts']
    Flux = '%.2e' %src['eflux']
    Flux_err = '%.2e' %src['eflux_err']
    Flux_UL = '%.2e' %src['eflux_ul95']
    Index = '%.2f' %(np.fabs(src['param_values'][1]))
    Index_err = '%.2f' %src['param_errors'][1]
    f = open('%s_Param.txt' % srcname, 'w')
    f.write(str(srcname) + "\t" + str(dist_sep) + "\t" + str(Flux) + "\t" + str(Flux_err) + "\t" + str(Index) + "\t" + str(Index_err) + "\t" + str(Flux_UL) + "\t" + str(TS) + "\n")
    f.close()
    
    return
    
################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":
    main(sys.argv)
