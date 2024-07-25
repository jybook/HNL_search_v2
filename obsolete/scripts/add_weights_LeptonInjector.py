#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1/icetray-start
#METAPROJECT /data/user/dvannerom/Metaprojects/OscNext/build

#==============================
# To use this script, one must source a /cvmfs environment containing LeptonWeighter,
# such as /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1/setup.sh
#==============================

import time
import os
from os.path import expandvars
import numpy as np
from I3Tray import *
from icecube import icetray, dataio, dataclasses, LeptonInjector
from icecube import phys_services, earthmodel_service
from icecube import tableio, hdfwriter
from optparse import OptionParser
from icecube.icetray import I3Units
import LeptonWeighter as LW

start_time = time.time()

usage = "usage: %prog [options] inputfile" 
parser = OptionParser(usage)
parser.add_option("-i","--infile", 		type="string", help = 'The infile Generation level file.',
                default="/data/user/dvannerom/HNL/test_Generation_level.i3.zst")
parser.add_option("-o","--outfile", 	type="string", help = 'The outfile name.',
                default="/data/user/dvannerom/HNL/test_Generation_level_weights.i3.zst")
parser.add_option("-c","--configuration",dest="configuration",type="string",
        #default = "/data/ana/SterileNeutrino/IC86/HighEnergy/HNL/MC/Tau/Gen/00001-01000/Generation_data.lic")
        default = "/data/ana/SterileNeutrino/IC86/HighEnergy/HNL/MC/M_1GeV/Gen/00001-01000/Generation_data.lic")
(options,args) = parser.parse_args() 

infile 		= options.infile
outfile     = options.outfile
#outfile 	 = '/data/user/dvannerom/HNL/TauGenWeightsFile/'+infile.split('/')[-1].replace('Gen','GenWeights')
#outfile 	= infile.replace('Gen','GenWeights')
#outfile     = outfile.replace('i3.zst','i3')
configuration = options.configuration

xs_location = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/xs_iso/'

tray = I3Tray()

#tray.AddModule("I3Reader", "reader",filenamelist=infile)
tray.Add("I3Reader",Filename=infile)

# Get generator-level weights
def get_weights(frame):
    #Flux = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'
    #Flux = '/home/carguelles/vault/scratch/DavidEmpire/script/data/conventional_atmospheric.hdf5'
    Flux = '/data/user/dvannerom/HNL/Fluxes/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5'
    LIConfiguration = configuration

    simulation_generation   = LW.MakeGeneratorsFromLICFile(LIConfiguration)
    nusquids_flux         = LW.nuSQUIDSAtmFlux(Flux)
    #dsdxdy_nu_CC = xs_location+'dsdxdy_nu_CC_iso.fits'
    #dsdxdy_nubar_CC = xs_location+'dsdxdy_nubar_CC_iso.fits'
    #dsdxdy_nu_NC = xs_location+'dsdxdy_nu_NC_iso.fits'
    #dsdxdy_nubar_NC = xs_location+'dsdxdy_nubar_NC_iso.fits'
    dsdxdy_nu_CC = "/home/carguelles/vault/scratch/DavidEmpire/script/data/dsdxdy-numu-N-cc-HERAPDF15NLO_EIG_central.fits"
    dsdxdy_nubar_CC = "/home/carguelles/vault/scratch/DavidEmpire/script/data/dsdxdy-numubar-N-cc-HERAPDF15NLO_EIG_central.fits"
    dsdxdy_nu_NC = "/home/carguelles/vault/scratch/DavidEmpire/script/data/dsdxdy-numu-N-nc-HERAPDF15NLO_EIG_central.fits"
    dsdxdy_nubar_NC = "/home/carguelles/vault/scratch/DavidEmpire/script/data/dsdxdy-numubar-N-nc-HERAPDF15NLO_EIG_central.fits"
    xs                      = LW.CrossSectionFromSpline(dsdxdy_nu_CC, dsdxdy_nubar_CC, dsdxdy_nu_NC, dsdxdy_nubar_NC)
    weighter            = LW.Weighter(nusquids_flux,xs,simulation_generation)

    LWevent = LW.Event()
    EventProperties                 = frame['EventProperties']
    LWevent.primary_type            = LW.ParticleType(EventProperties.initialType)
    LWevent.final_state_particle_0  = LW.ParticleType(EventProperties.finalType1)
    LWevent.final_state_particle_1  = LW.ParticleType(EventProperties.finalType2)
    LWevent.zenith                  = EventProperties.zenith
    LWevent.energy                  = EventProperties.totalEnergy
    LWevent.azimuth                 = EventProperties.azimuth
    LWevent.interaction_x           = EventProperties.finalStateX
    LWevent.interaction_y           = EventProperties.finalStateY
    LWevent.total_column_depth      = EventProperties.totalColumnDepth
    LWevent.radius                  = EventProperties.impactParameter
    LWevent.x = 0.
    LWevent.y = 0.
    LWevent.z = 0.

    weight = weighter(LWevent)/2.
    frame.Put("weight",dataclasses.I3Double(weight))

#tray.AddModule(get_weights,streams=[icetray.I3Frame.DAQ])
#tray.AddModule("I3NullSplitter","nullsplit")
tray.Add(get_weights,streams=[icetray.I3Frame.DAQ])

tray.Add("I3Writer",Filename=outfile)

tray.Execute()
tray.Finish()

print('done ...')
