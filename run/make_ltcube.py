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


def main(cmd_line):

    '''Here we can set variable for path definitions
       In addition to all desired parameters, 
       double check you have the correct LTcube, extdir
    '''

####################################
	#Get inputs
####################################S
    srcname = cmd_line[1]
    ra      = float(cmd_line[2])
    dec     = float(cmd_line[3])

####################################
    #Create config file
####################################
    ft1 = '/Users/asteinhe/FermiLAT/BHinEGs_DM/data/'+srcname+'/ft1.fits' #data files
    ft2 = '/Users/asteinhe/FermiLAT/BHinEGs_DM/data/ft2_sc.fits' #SC file

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
        yml.write("#--------#\n")
        yml.write("binning:\n")
        yml.write("  roiwidth : 8\n")
        yml.write("  binsz : 0.08\n")
        yml.write("  binsperdec : 4\n")
        yml.write("#--------#\n")
        yml.write("selection:\n")
        yml.write("  emin : 500\n")
        yml.write("  emax : 1000000\n")
        yml.write("  zmax : 100\n")
        yml.write("  target : '%s'\n" % srcname)
        yml.write("  radius : 15\n")
        yml.write("  tmin : 239557417\n")
        yml.write("  tmax : 665410983\n")
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
        yml.write("  src_radius : 10\n")
        yml.write("  src_roiwidth : 10\n")
        yml.write("  galdiff : '%s'\n" %galdiff)
        yml.write("  isodiff : '%s'\n" %isodiff)
        yml.write("  catalogs :\n")
        yml.write("    - '4FGL-DR3'\n")
        yml.write("  sources :\n")
        yml.write("    - { 'name' : '%s', 'ra' : %s, 'dec' : %s, 'SpectrumType' : PowerLaw }\n" % (srcname,ra,dec))
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
    
    return
    
################################################################################################
################################################################################################
################################################################################################
if __name__=="__main__":
    main(sys.argv)
