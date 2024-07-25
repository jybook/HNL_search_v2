#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Use with OscNext meta-project (e.g. Leander's):
# unset OS_ARCH; eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh
# export PYTHONPATH=$PYTHONPATH:/afs/ifh.de/group/amanda/scratch/lfischer/software/fridge
# /afs/ifh.de/user/l/lfischer/scratch/ICsoft/meta-projects/oscnext_meta/build/env-shell.sh

# OscNext files:
# gcd: /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz

# nominal sets L5 at zeuthen:
# genie level5_v01.04 (a.trettin): /lustre/fs23/group/icecube/trettin/oscNext/processing/pass2/genie/level5_v01.04

# nominal sets L7 at zeuthen:
# genie level7_v02.00 : /lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/genie/level7_v02.00/1[2,4,6]0000
# muongun level7_v02.00 : /lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/muongun/level7_v02.00/130000
# noise level7_v02.00 : /lustre/fs22/group/icecube/lfischer/data/OscNext/pass2/noise/level7_v02.00/888003

# HNL signal files at zeuthen:
# MC (different sets): /lustre/fs22/group/icecube/lfischer/data/HNL/SterileNeutrino/IC86/HighEnergy/HNL/MC

#
# Author: Leander Fischer
#

import time
t0 = time.time()

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-i", "--infile",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/SterileNeutrino/IC86/HighEnergy/HNL/MC/M_1GeV/Ares/IC86.AVG/L7/domeff_0.97/00001-01000/L7_00_11_00056.i3.zst",
                 dest="INFILE", help="Read input from INFILE (.i3{.gz/.zst} format)")
parser.add_option("-o", "--outfile",
                #  default="/lustre/fs22/group/icecube/lfischer/data/HNL/test_output_all.i3.zst",
                 dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz/.zst} format)")
parser.add_option("-g", "--gcdfile", default="/lustre/fs22/group/icecube/lfischer/data/HNL/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
                 dest="GCDFILE", help="Read in GCD file")
parser.add_option("-v", "--verbose", default=False, action="store_true",
                 dest="VERBOSE", help="print info level output")
parser.add_option("-l", "--logtofile", default=False, action="store_true",
                 dest="LOGTOGILE", help="Store log into a file with the current datetime string as name.")
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

if options.LOGTOGILE:
    timestring = time.strftime("%Y%m%d_%H%M%S_monopod_and_taupede.log")
    logpath = os.path.join('/lustre/fs22/group/icecube/lfischer/LOGS/monopod_and_taupede', timestring)
    logging.rotating_files(logpath)

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


# # Check if the event was contained in deepcore.
# def DCContainment(frame):
#     '''
#     Define flag indicating if the event (DIS vertex) is contained in DC volumes
#     '''
#     logging.log_info('Checking for DC containment.')

#     # Get the primary
#     mctree = frame['I3MCTree']
#     true_particle = mctree.primaries[0]
#     true_vertex = true_particle.pos
    
#     # DeepCore containment
#     deepcore_containment_bool = bool(np.sqrt( np.square(true_vertex.x - 56.290) + np.square(true_vertex.y - -34.880) ) < 145.) # Radial distance from string 36
#     deepcore_containment_bool = bool(deepcore_containment_bool and (true_vertex.z <= -201.93 ) and (true_vertex.z >= -476.93 )) # Vertical region
#     # Add the flags to the frame

#     frame.Put("DCContainment", I3Bool(deepcore_containment_bool))
#     return True


# # Calculate true decay lenght and add it to the frame (for quick testing).
# def TrueDecayLength(frame):
#     '''
#     Calculate true decay lenght from MC hadron objecs.
#     '''
#     logging.log_info('Calculating true decay length.')

#     mctree = frame['I3MCTree']        
#     p_true = mctree.primaries[0]

#     p_daughters = mctree.get_daughters(p_true)

#     assert(len(p_daughters) == 2)

#     for p_daughter in p_daughters:
#         if p_daughter.type == dataclasses.I3Particle.Hadrons:
#             casc_0_true = p_daughter
#     # MODIFICATION (for 190605 set, which has two hadron objects and intermediate HNL particle)
#         else:
#             casc_1_true = mctree.get_daughters(p_daughter)[0]

#     decay_l_true = (phys_services.I3Calculator.distance(casc_0_true,casc_1_true)/icetray.I3Units.m)

#     frame.Put("TrueDecayLength", I3Double(decay_l_true))
#     return True


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


#######
###   TAUPEDE
########
@icetray.traysegment
def TaupedeFitWrapper(tray, name, Seed, Iterations, PhotonsPerBin,
                      ethreshold, outerboundary_xy, outerboundaries_z,
                      lengthbounds, **millipede_params):
    
    logging.log_info('Running Taupede fit.')


    ############ TAUPEDE VERSIONS ############

    # full fit v0     (timed) full fit of all parameters in one go
    # split fit v2    first fitting length and energies then full fit
    # split fit v3    first fitting length+energies+time and then full fit
    # split fit v4    first untimed fit of all parametes but the time and then full fit

    # split fit v2 minuit2 simplex      same with minuit2 simplex
    # split fit v2 minuit2 migrad       same with minuit2 migrad
    # split fit v2 length seeds         same with five different length seeds

    """
    single seed taupede fit (original fit chain )
    run:no
    """

    # amptag = name+'AmpFit'
    # recotag = name

    # # amp fit    
    # tray.Add(TaupedeFit,
    #         amptag,
    #         Seed = Seed,
    #         PhotonsPerBin = -1,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         StepT = 0,
    #         Iterations = Iterations,
    #         **millipede_params)
    
    # # timed fit (full)    
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         # Seed = amptag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    single seed llh scan fit chain
    run:no
    """

    # recotag = name

    # # llh scan fit chain
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepT=0,
    #         StepD=0,
    #         StepZenith=0,
    #         StepAzimuth=0,
    #         # StepL=0,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)


    """
    full fit v0 (timed fit (full))
    run:yes
    """

    # recotag = name
   
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    full fit v1 (timed fit (full) 4x iterations)
    run:no
    """

    # recotag = name
   
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = 4,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    split fit v0 (first fit only length (and energies) then fit full)  PHOTONSPERBIN was not set to -1 (check if v2 is the same)
    run:no
    """

    # recotag = name
    # length_fit_tag = name+'LengthFit'

    # # length fit (not timed)
    # tray.Add(TaupedeFit, 
    #         length_fit_tag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepT=0,
    #         StepD=0,
    #         StepZenith=0,
    #         StepAzimuth=0,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    # # full fit
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = length_fit_tag,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    split fit v2 (same as split fit v0, but with PhotonsPerBin=-1 (to set it to really non timed)) (first fit only length (and energies) then fit full)
    run:yes
    """

    recotag = name
    length_fit_tag = name+'LengthFit'

    # length fit (not timed)
    tray.Add(TaupedeFit, 
            length_fit_tag,
            Seed = Seed,
            LengthBounds = lengthBounds,
            StepT=0,                    # untimed!
            StepD=0,
            StepZenith=0,
            StepAzimuth=0,
            StepL = 1,
            Iterations = Iterations,
            PhotonsPerBin = -1,          # untimed!
            **millipede_params)

    # full fit
    tray.Add(TaupedeFit, 
            recotag,
            Seed = length_fit_tag,
            LengthBounds = lengthBounds,
            StepL = 1,
            Iterations = Iterations,
            PhotonsPerBin = PhotonsPerBin,
            **millipede_params)

    """
    split fit v3 (same as split fit v2, but this time running a length+time fit first) (first fit length+time (and energies) then fit full)
    run:yes
    """

    # recotag = name
    # length_fit_tag = name+'LengthFit'

    # # length fit (not timed)
    # tray.Add(TaupedeFit, 
    #         length_fit_tag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepD=0,
    #         StepZenith=0,
    #         StepAzimuth=0,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    # # full fit
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = length_fit_tag,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    split fit v4 (original fit chain with 1x iterations)
    run:yes
    """

    # amptag = name+'AmpFit'
    # recotag = name

    # # amp fit    
    # tray.Add(TaupedeFit,
    #         amptag,
    #         Seed = Seed,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         StepT = 0,
    #         Iterations = Iterations,
    #         PhotonsPerBin = -1,
    #         **millipede_params)
    
    # # timed fit (full)    
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = amptag,
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    split fit v2 length seeds (same as split fit v2, but brute force length seed for length fit and then select best)
    run:yes
    """

    # scales = [5., 25., 50., 100., 200.]

    # def addbruteforceseeds(frame, Seed, Output):

    #     logging.log_info('Adding brute force seeds to the frame.')

    #     taupedeseedbase = I3Particle(frame[Seed])
    #     # print(taupedeseedbase)
    #     for scale in scales:
    #         taupedeseed = I3Particle(taupedeseedbase)
    #         taupedeseed.length = scale
    #         # print(taupedeseed)
    #         frame.Put(Output+'_scale%i' % scale, taupedeseed)
    # tray.Add(addbruteforceseeds, 'AddBruteForceSeeds', Seed=Seed, Output='TaupedeBruteForceWrapperSeed')

    # for scale in scales:
    #     seedtag = 'TaupedeBruteForceWrapperSeed_scale%i' % (scale)
    #     length_fit_tag = 'TaupedeBruteForceWrapperReco_scale%i' % (scale)

    #     # length fit (not timed)
    #     tray.Add(TaupedeFit, 
    #             length_fit_tag,
    #             Seed = seedtag,
    #             LengthBounds = lengthBounds,
    #             StepT=0,                    # untimed!
    #             StepD=0,
    #             StepZenith=0,
    #             StepAzimuth=0,
    #             StepL = 1,
    #             Iterations = Iterations,
    #             PhotonsPerBin = -1,          # untimed!
    #             **millipede_params)
        
        
    # def findbestfit(frame, Output):

    #     logging.log_info('Selecting best length fit to seed the full fit.')

    #     rlogl_map = []

    #     for scale in scales:
    #         length_fit_tag = 'TaupedeBruteForceWrapperReco_scale%i' % (scale)

    #         fitparamstag = length_fit_tag+'FitParams'

    #         if not fitparamstag in frame:
    #             continue

    #         rlogl = frame[fitparamstag].rlogl

    #         rlogl_map.append((length_fit_tag, rlogl))

    #     # print rlogl_map

    #     dtype = [('fittag', 'S60'), ('rlogl', float)]
    #     fitmap = np.array(rlogl_map, dtype=dtype)

    #     sortedfitmap = np.sort(fitmap, order='rlogl')
    
    #     bestfittag = sortedfitmap[0]['fittag']

    #     bestfitparticle = I3Particle(frame[bestfittag])

    #     # print(bestfitparticle)

    #     frame.Put(Output, bestfitparticle) 

    # tray.Add(findbestfit, 'FindBestFit', Output='TaupedeBestLengthFit')

    # recotag = name
  
    # # full fit
    # tray.Add(TaupedeFit, 
    #         recotag,
    #         Seed = 'TaupedeBestLengthFit',
    #         LengthBounds = lengthBounds,
    #         StepL = 1,
    #         Iterations = Iterations,
    #         PhotonsPerBin = PhotonsPerBin,
    #         **millipede_params)

    """
    determine a good fit using the particle properties instead of the rlogl
    """
    def isgoodfit(cascade0, cascade1):
        logging.log_info('Checking for fit status(geometry/energy).')
        goodfit = True
        # check for reco energies (E0, E1) are above ethreshold
        if cascade0.energy < ethreshold or cascade1.energy < ethreshold:
            goodfit = False
        # check for containment of cascade 0
        if np.abs(cascade0.pos.x) > outerboundary_xy or np.abs(cascade0.pos.y) > outerboundary_xy or outerboundaries_z[0] > cascade0.pos.z or cascade0.pos.z > outerboundaries_z[1]:
            goodfit = False
        # check for containment of cascade 1
        if np.abs(cascade1.pos.x) > outerboundary_xy or np.abs(cascade1.pos.y) > outerboundary_xy or outerboundaries_z[0] > cascade1.pos.z or cascade1.pos.z > outerboundaries_z[1]:
            goodfit = False
        return goodfit


    """
    add manual fit status information of the final best fit
    """
    def addfitstatus(frame, name, Seed):
        logging.log_info('Adding seed fit status.')
        fitstatus = I3Particle.FitStatus.OK
        outtag = name + 'ManualFitStatus'

        # check for Seed Particle Seed
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

# pulses='SRTTWOfflinePulsesDC'
pulses='SplitInIcePulses'

tray = I3Tray()

tray.context['I3FileStager'] = dataio.get_stagers()

tray.AddModule('I3Reader',
               'reader',
               FilenameList=[options.GCDFILE, options.INFILE],
              )

# tray.AddModule(DCContainment, "DCContainment", Streams=[icetray.I3Frame.Physics])

# tray.AddModule(TrueDecayLength, "TrueDecayLength", Streams=[icetray.I3Frame.Physics])

tray.AddModule(checky, keys=[pulses])

# set taupede Seed particle from retro
tray.AddModule(CreateRetroSeed, 'CreateRetroSeed', Streams=[icetray.I3Frame.Physics])

# # add table base/cascade service (old)
# table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/ems_spice1_z20_a10.%s.fits')
# cascade_service = photonics_service.I3PhotoSplineService(table_base % 'abs', table_base % 'prob', 0., maxRadius=480.)

# add table base/cascade service (newest photonics tables including effective distande parametrization by Marcel Usner)
# needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
# splinetabledir = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines'
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
    # minimizer='MIGRAD'
    # minimizer='LBFGSB'
    )

# logging.log_info('taupede_params: {}'.format(taupede_params))

PhotonsPerBin = 1
ethreshold = 0.1
outerboundary_xy = 250.
outerboundaries_z = [-600., 100.]
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
        outerboundary_xy = outerboundary_xy,
        outerboundaries_z = outerboundaries_z,
        lengthBounds = lengthBounds,
        **taupede_params)

tray.Add("I3Writer", Streams=map(icetray.I3Frame.Stream, 'QP'),
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=options.OUTFILE)

tray.Execute()
tray.Finish()

t1 = time.time()
logging.log_info('Time it took: {:.3f} s'.format(t1-t0))