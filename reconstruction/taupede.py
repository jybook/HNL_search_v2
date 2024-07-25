#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
(options,args) = parser.parse_args()
if not options.INFILE:
    parser.error('Infile not specified.')
if not options.GCDFILE:
    parser.error('GCDfile not specified.')


import os
import numpy as np
from icecube import icetray, dataclasses, dataio
from I3Tray import *
from icecube.dataclasses import I3Double, I3Particle, I3Constants
from icecube.icetray import I3Bool

from icecube import photonics_service

from icecube.millipede import TaupedeFit

import icecube.icetray.i3logging as logging


if options.VERBOSE:
    logging.set_level('INFO')
else:
    logging.set_level('WARN')

logging.log_info('Options: {}'.format(options))

OUTPATH = os.path.split(options.OUTFILE)[0]
if not os.path.exists(OUTPATH):
    logging.log_info('Creating output filepath: {}'.format(OUTPATH))
    os.mkdir(OUTPATH)

# Log length of input pulsemap.
def checky(frame, keys):
    for k in keys:
        # print(k, len(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, k)))
        logging.log_info('Input Pulses (length) : {} ({})'.format(k, len(dataclasses.I3RecoPulseSeriesMap.from_frame(frame, k))) )


# Create seed particle from Retro best fit results.
def CreateRetroSeed(frame, seedname='RetroSeedParticle'):
    '''
    Create seed particle from retro best fit particle with em cascade energy
    '''
    logging.log_info('Creating retro reco seed.')

    retro_key_base = 'L7_reconstructed_'

    seed_particle = I3Particle()
    # set vertex
    seed_particle.pos = dataclasses.I3Position(
        frame[retro_key_base+'vertex_x'].value,
        frame[retro_key_base+'vertex_y'].value,
        frame[retro_key_base+'vertex_z'].value,
    )
    
    # set direction
    seed_particle.dir = dataclasses.I3Direction(
        frame[retro_key_base+'zenith'].value,
        frame[retro_key_base+'azimuth'].value,
    )

    # set track length
    seed_particle.length = frame[retro_key_base+'track_length'].value
    
    # set energy as em cascade energy
    seed_particle.energy = frame[retro_key_base+'em_cascade_energy'].value

    # set time
    seed_particle.time = frame[retro_key_base+'time'].value

    # set fit status
    seed_particle.shape = I3Particle.Primary
    seed_particle.fit_status = I3Particle.OK

    # set particle type
    seed_particle.type = I3Particle.NuTau
    seed_particle.pdg_encoding = 16

    seed_particle.speed = dataclasses.I3Constants.c

    # print(seed_particle)

    frame.Put(seedname, seed_particle)
    return True


###########################
###   TAUPEDE
###########################

@icetray.traysegment
def TaupedeFitWrapper(tray, name, Seed, Iterations, PhotonsPerBin,
                      ethreshold, radial_boundary, vertical_outerboundaries,
                      lengthbounds, **millipede_params):
    
    logging.log_info('Running Taupede fit.')

    """
    Chosen fit routine: Fit is split into length+energies fit and full fit,
    where the latter is using the result from the first as a seed.
    Length+energies fit are additionally seeded with 3 scaled values of the
    Retro track length as seed. Best performing is selected.
    """

    scales = [0.5, 1., 1.5]

    def addbruteforceseeds(frame, Seed, Output):

        logging.log_info('Adding brute force seeds to the frame.')

        taupedeseedbase = I3Particle(frame[Seed])
        # print(taupedeseedbase)
        for scale in scales:
            taupedeseed = I3Particle(taupedeseedbase)
            taupedeseed.length = scale * taupedeseed.length
            # print(taupedeseed)
            frame.Put(Output+'_scale{:.01f}'.format(scale), taupedeseed)

    tray.Add(addbruteforceseeds, 'AddBruteForceSeeds', Seed=Seed, Output='TaupedeBruteForceWrapperSeed')

    for scale in scales:
        seedtag = 'TaupedeBruteForceWrapperSeed_scale{:.01f}'.format(scale)
        length_fit_tag = 'TaupedeBruteForceWrapperReco_scale{:.01f}'.format(scale)

        # length fit (not timed)
        tray.Add(TaupedeFit, 
                length_fit_tag,
                Seed = seedtag,
                LengthBounds = lengthBounds,
                StepT=0,
                StepD=0,
                StepZenith=0,
                StepAzimuth=0,
                StepL = 1,
                Iterations = Iterations,
                PhotonsPerBin = -1,
                **millipede_params)


    def findbestfit(frame, Seed, Output):

        logging.log_info('Selecting best length fit to seed the full fit.')

        rlogl_map = []

        for scale in scales:
            length_fit_tag = 'TaupedeBruteForceWrapperReco_scale{:.01f}'.format(scale)

            fitparamstag = length_fit_tag+'FitParams'

            if not fitparamstag in frame:
                continue

            rlogl = frame[fitparamstag].rlogl

            rlogl_map.append((length_fit_tag, rlogl))

        # print rlogl_map

        dtype = [('fittag', 'S60'), ('rlogl', float)]
        fitmap = np.array(rlogl_map, dtype=dtype)

        sortedfitmap = np.sort(fitmap, order='rlogl')
    
        if not len(sortedfitmap) == 0:
            bestfittag = sortedfitmap[0]['fittag']

            bestfitparticle = I3Particle(frame[bestfittag])
        else:
            bestfitparticle = I3Particle(frame[Seed])

        # print(bestfitparticle)

        frame.Put(Output, bestfitparticle) 

    tray.Add(findbestfit, 'FindBestFit', Seed=Seed, Output='TaupedeBestLengthFit')

    recotag = name
  
    # full fit
    tray.Add(TaupedeFit, 
            recotag,
            Seed = 'TaupedeBestLengthFit',
            LengthBounds = lengthBounds,
            StepL = 1,
            Iterations = Iterations,
            PhotonsPerBin = PhotonsPerBin,
            **millipede_params)


    """
    determine a good fit using the particle properties energy and soft containment
    """
    def isgoodfit(cascade0, cascade1):
        logging.log_info('Checking for fit status(geometry/energy).')
        goodfit = True
        # check for reco energies (E0, E1) are above ethreshold
        if cascade0.energy < ethreshold or cascade1.energy < ethreshold:
            goodfit = False
        # check for containment of cascade 0
        cascade0_radial = (np.sqrt( np.square(cascade0.pos.x - 56.290) + np.square(cascade0.pos.y - -34.880) ) > radial_boundary) # distance to string 36
        cascade0_depth = (cascade0.pos.z >= vertical_outerboundaries[1]) + (cascade0.pos.z <= vertical_outerboundaries[0] ) # Vertical region
        if cascade0_radial or cascade0_depth:
            goodfit = False
        # check for containment of cascade 1
        cascade1_radial = (np.sqrt( np.square(cascade1.pos.x - 56.290) + np.square(cascade1.pos.y - -34.880) ) > radial_boundary) # distance to string 36
        cascade1_depth = (cascade1.pos.z >= vertical_outerboundaries[1]) + (cascade1.pos.z <= vertical_outerboundaries[0] ) # Vertical region
        if cascade1_radial or cascade1_depth:
            goodfit = False
        return goodfit


    """
    add manual fit status information of the final best fit
    """
    def addfitstatus(frame, name, Seed):
        logging.log_info('Adding seed fit status.')
        fitstatus = I3Particle.FitStatus.OK
        outtag = name + 'ManualFitStatus'

        # check if Particle "Seed" exists
        if Seed not in frame:
            fitstatus = I3Particle.FitStatus.MissingSeed
            frame.Put(outtag, I3Double(fitstatus))
            return True
        
        # check fit result in name/nameParticles
        if name not in frame or name + 'Particles' not in frame:
            fitstatus = I3Particle.FitStatus.GeneralFailure
            frame.Put(outtag, I3Double(fitstatus))
            return True

        # check for general weirdness in the fit result        
        cascade0, cascade1 = frame[name + 'Particles']        
        if (cascade0.energy + cascade1.energy == 0.) or cascade0.pos.r < 0. or cascade1.pos.r < 0.:
            fitstatus = I3Particle.FitStatus.FailedToConverge
            frame.Put(outtag, I3Double(fitstatus))
            return True
            
        # check for soft containment and energy threshold
        goodfit = isgoodfit(cascade0, cascade1)
        if not goodfit:
            fitstatus = I3Particle.FitStatus.InsufficientQuality
            frame.Put(outtag, I3Double(fitstatus))
            return True

        frame.Put(outtag, I3Double(fitstatus))
        return True

    tray.Add(addfitstatus, name=name, Seed=Seed)

###########################


"""
Perform reconstruction and pre-processing
"""

# use uncleaned pulses since millipede accounts for noise
pulses='SplitInIcePulses'

tray = I3Tray()
    
tray.context['I3FileStager'] = dataio.get_stagers()

tray.AddModule('I3Reader',
               'reader',
               FilenameList=[options.GCDFILE, options.INFILE],
              )

tray.AddModule(checky, keys=[pulses])

# set taupede seed particle from retro
tray.AddModule(CreateRetroSeed, 'CreateRetroSeed', Streams=[icetray.I3Frame.Physics])

# add table base/cascade service (newest photonics tables including effective distande parametrization by Marcel Usner)
# needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
                                                         timingtable = table_base.format('prob'),
                                                         effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
                                                         timingSigma = 0,
                                                         maxRadius=400.
                                                        )

# set TaupedeFit parameters
taupede_params = dict(
    Pulses=pulses,
    CascadePhotonicsService=cascade_service,
    minimizer='SIMPLEX'
    )

# logging.log_info('taupede_params: {}'.format(taupede_params))

PhotonsPerBin = 1
ethreshold = 0.0
# DC spatial boundareis
radial_boundary = 150.
vertical_outerboundaries = [-500., -150.]
lengthBounds = [-800., 800.]

name = 'TaupedeFit'
Seed = 'RetroSeedParticle'

# run the fit (through the TaupedeFitWrapper)
tray.Add(TaupedeFitWrapper,
        name,
        Seed = Seed,
        Iterations = 1,
        PhotonsPerBin = PhotonsPerBin,
        ethreshold = ethreshold,
        radial_boundary = radial_boundary,
        vertical_outerboundaries = vertical_outerboundaries,
        lengthBounds = lengthBounds,
        **taupede_params)

tray.Add("I3Writer",
    # Streams=map(icetray.I3Frame.Stream, 'QP'),
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=options.OUTFILE)

tray.Execute()
tray.Finish()

t1 = time.time()
logging.log_info('Time it took: {:.3f} s'.format(t1-t0))