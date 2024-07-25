#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/


### The MuMillipede reconstruction script ###

#
# Author: Leander Fischer
#

import time
t0 = time.time()

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--infile",
                 dest="INFILE", help="Read input from INFILE (.i3{.gz/.zst} format)")
parser.add_option("-o", "--outfile",
                 dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz/.zst} format)")
parser.add_option("-g", "--gcdfile", default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
                 dest="GCDFILE", help="Read in GCD file")
parser.add_option("-p", "--pulseseries", default='SplitInIcePulses',
                 dest="PULSESERIES", help="Choose pulse series to use")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
parser.add_option("--hdf5", action="store_false",
                 dest="write_hdf5", default=True, help="Also store hdf5 file")
(options,args) = parser.parse_args()
if not options.INFILE:
    parser.error('Infile not specified.')
if not options.GCDFILE:
    parser.error('GCDfile not specified.')


write_hdf5  = options.write_hdf5


import os
import numpy as np
from icecube import icetray, dataclasses, dataio, photonics_service, phys_services
from I3Tray import *
from icecube.dataclasses import I3Double, I3Particle, I3Constants
from icecube.icetray import I3Bool

from icecube.millipede import TaupedeFit

import icecube.icetray.i3logging as logging

from copy import copy
from icecube.hdfwriter import I3HDFWriter


if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

logging.log_info('Options: {}'.format(options))

OUTPATH = os.path.split(options.OUTFILE)[0]
if not os.path.exists(OUTPATH):
    logging.log_info('Creating output filepath: {}'.format(OUTPATH))
    os.makedirs(OUTPATH)


def flip_cascades(cascades):
    '''
    Function to flip the physical properties of the taupede 
    cascade objects and also the direction/decay lengths.
    '''    
#     flip cascade directions, lengths
    for cascade in cascades:
        if cascade.length < 0.0:cascade.length = cascade.length * (-1)

#     exchange cascade positions
    tmp_pos_casc1 = copy(cascades[1].pos)
    cascades[1].pos = cascades[0].pos
    cascades[0].pos = tmp_pos_casc1

#     exchange cascade times
    tmp_time_casc1 = copy(cascades[1].time)
    cascades[1].time = cascades[0].time
    cascades[0].time = tmp_time_casc1

#     exchange cascade energies
    tmp_time_cenergy = copy(cascades[1].energy)
    cascades[1].energy = cascades[0].energy
    cascades[0].energy = tmp_time_cenergy

def flip_lengths_position_time(primary, cascades):
    '''
    Function to flip direction/decay length of a primary particle.
    '''    
    primary.length = primary.length * (-1)
    
    primary.pos = cascades[0].pos
    primary.time = cascades[0].time


def correct_negative_decay_lengths(frame, fitname='TaupedeFit'):
    '''
    Function to correct for negative reconstructed decay lenghts/flipped cascade order.
    '''    
    fit_particle_key = fitname
    fit_cascades_key = fit_particle_key + 'Particles'

    if frame.Has(fit_cascades_key):
        primary = frame[fit_particle_key]
        cascades = frame[fit_cascades_key]
        
        if(primary.length < 0.0):        
            flip_cascades(cascades)
            flip_lengths_position_time(primary, cascades)

            frame.Delete(fit_particle_key)
            frame.Put(fit_particle_key, primary)

            frame.Delete(fit_cascades_key)
            frame.Put(fit_cascades_key, cascades)
        
        # add reco cascades, energy and decay length to the frame
        frame['reco_decay_length'] = dataclasses.I3Double(primary.length)
        frame['reco_total_energy'] = dataclasses.I3Double(primary.energy)
        frame['reco_cascade_0'] = cascades[0]
        frame['reco_cascade_1'] = cascades[1]

    return True

# def write_true_data_to_keys(frame):
#     '''
#     Function to write true cascade information, energy and decay lenght to frame
#     '''    
#     cascades = frame['I3MCTree'].get_daughters(frame['I3MCTree'].primaries[0])

#     # add true cascades, energy and decay length to the frame
#     frame['true_decay_length'] = dataclasses.I3Double(phys_services.I3Calculator.distance(cascades[0], cascades[1]))
#     frame['true_total_energy'] = dataclasses.I3Double(cascades[0].energy + cascades[1].energy)
#     frame['true_cascade_0'] = cascades[0]
#     frame['true_cascade_1'] = cascades[1]

#     return True


"""
Perform reconstruction and pre-processing
"""

# use uncleaned pulses since millipede accounts for noise
pulses = options.PULSESERIES

tray = I3Tray()
    
tray.context['I3FileStager'] = dataio.get_stagers()

tray.AddModule('I3Reader',
               'reader',
               FilenameList=[options.GCDFILE, options.INFILE],
              )

# add table base/cascade service (photonics tables including effective distande parametrization by Marcel Usner)
# needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
                                                         timingtable = table_base.format('prob'),
                                                         effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
                                                         timingSigma = 0,
                                                         maxRadius=400.
                                                        )

# # try with the new photonics tables that tianlu produced (end of 2021)
# # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
# table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_3.2.1_flat_z20_a5.{}.fits')
# cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
#                                                          timingtable = table_base.format('prob'),
#                                                          effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_3.2.1_z20.eff.fits'),
#                                                          tiltTableDir = os.path.expandvars( '$I3_SRC/photonics-service/resources/tilt/'),
#                                                          timingSigma = 0,
#                                                          maxRadius=400.
#                                                         )

# set MuMillipedeFit parameters
common_params = dict(
    Pulses=pulses,
    CascadePhotonicsService=cascade_service,
    # minimizer='SIMPLEX',
    PhotonsPerBin = 1,
    )

# logging.log_info('common_params: {}'.format(common_params))

name = 'MuMillipedeFit'
Seed = 'TaupedeFit'

# run function to correct for negative decay lenghts and write reco data as keys
tray.AddModule(correct_negative_decay_lengths, 'CorrectNegativeDecayLength', fitname=Seed)

# # add function to write true data as keys
# tray.AddModule(write_true_data_to_keys, 'WriteTrueDataToKeys')

# run the fit (with the TaupedeFit as Seed)
tray.AddModule("MuMillipede",
        name,
        Output=name,
        ShowerSpacing=5,
        MuonSpacing=0,
        SeedTrack = Seed,
        **common_params)


tray.Add("I3Writer",
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=options.OUTFILE)


if write_hdf5:

    hdf5_keys = [
        'true_decay_length',
        'true_total_energy',
        'true_cascade_0',
        'true_cascade_1',
        'reco_decay_length',
        'reco_total_energy',
        'reco_cascade_0',
        'reco_cascade_1',
        'PoleMuonLlhFit',
        'TaupedeFit',
        'TaupedeFit'+'FitParams',
        # 'TaupedeFit'+'ManualFitStatus',
    ]

    hdf5_keys.append(name)

    logging.log_info('hdf5 keys to store: {}'.format(hdf5_keys))

    outfile_hdf5 = options.OUTFILE.replace('.i3.zst',".hdf5")

    tray.Add(I3HDFWriter,
        Output=outfile_hdf5,
        Keys=hdf5_keys,
        SubEventStreams=["InIceSplit"]
    )


tray.Execute()
tray.Finish()

t1 = time.time()
logging.log_info('Time it took: {:.3f} s'.format(t1-t0))
