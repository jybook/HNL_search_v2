#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description:
This python script is used to make 1-d, 2-d likelihood surface plots. It uses the gulliview.py module from gulliver_modules, which i modified for my purposes. In this case i selected a specific event i isolated beforehand.

Plot likelihood space around true mc parameters (probably not the minimum but we will see).
"""

import os

import numpy as np

import icecube
from icecube import icetray, dataclasses, dataio, phys_services, photonics_service

from icecube.dataclasses import I3Double, I3Particle, I3VectorI3Particle, I3Constants

import icecube.gulliver_modules
import icecube.gulliver_modules.gulliview

import icecube.lilliput
import icecube.lilliput.segments
import icecube.linefit

from icecube.millipede import MonopodFit, TaupedeFit, MillipedeFitParams, HighEnergyExclusions

from icecube.icetray import I3Units
import icecube.icetray.i3logging as logging

from I3Tray import I3Tray

import time

t0 = time.time()

def SelectEvent(frame):
    '''
    Define flag that only selects single event with given energy/zenith
    '''

    logging.log_info('Searching wanted event.')

    # Get the primary
    mctree = frame['I3MCTree']
    true_particle = mctree.primaries[0]
    
    # # this event is in file L7_00_11_00056.i3.zst
    # chosen_energy = 75.83558806963663
    # chosen_zenith = 2.53893805725283

    # # event 1 in file L7_00_11_01191.i3.zst of 190605
    # chosen_energy = 99.03545958187648
    # chosen_zenith = 3.0184764924536753

    # # event 2 in file L7_00_11_01791.i3.zst of 190605
    # chosen_energy = 133.45304294720898
    # chosen_zenith = 3.0015342375731393    

    # event 4 in file L7_00_11_02153.i3.zst of 190605
    chosen_energy = 30.894653597352505
    chosen_zenith = 2.9970802134791934

    selection_bool = (true_particle.energy == chosen_energy) and (true_particle.dir.zenith == chosen_zenith)
    if selection_bool:
        logging.log_info('Found it, here we go.')

    return bool(selection_bool)


# def ContainmentFlags(frame):
#     '''
#     Define flags indicating if the event is contained in fiducial volumes
#     '''
#     # Get the primary
#     mctree = frame['I3MCTree']
#     true_particle = mctree.primaries[0]
#     # Check this particle is a neutrino, othewise containment does not apply
#     if np.abs(true_particle.pdg_encoding) not in [12,14,16]:
#         return False
#     # Upgrade containment (v53 geom)
#     true_vertex = true_particle.pos
#     # upgrade_containment_bool = np.sqrt( np.square(true_vertex.x - 47.29) + np.square(true_vertex.y - -57.0125) ) <= 50. # Radial distance from string 91, v53 geom
#     # upgrade_containment_bool = upgrade_containment_bool and (true_vertex.z <= -201.93 ) and (true_vertex.z >= -476.93 )
#     # DeepCore containment
#     deepcore_containment_bool = np.sqrt( np.square(true_vertex.x - 56.290) + np.square(true_vertex.y - -34.880) ) < 145. # Radial distance from string 36
#     deepcore_containment_bool = deepcore_containment_bool and (true_vertex.z <= -201.93 ) and (true_vertex.z >= -476.93 ) # Vertical region
#     # Add the flags to the frame
#     containment_map = dataclasses.I3MapStringBool()
#     # containment_map["Upgrade"] = bool(upgrade_containment_bool)
#     containment_map["DeepCore"] = bool(deepcore_containment_bool)
#     frame.Put("ContainmentFlags",containment_map)
            
#     return bool(deepcore_containment_bool)   


class StoreMCTruth(icecube.icetray.I3ConditionalModule):
    """
    Simple module to store the MC Truth information to a text file.
    Set up as I3ConditionalModule so it can iterate over the frame numbers.
    """

    @property
    def Filename(self):
        """
        Output filename base for plots; if None, run interactively.
        """
        pass

    def __init__(self, ctx):
        icecube.icetray.I3ConditionalModule.__init__(self, ctx)
        self.filename = None
        self.AddParameter("Filename",
                          self.__class__.Filename.__doc__, self.filename)

    def Configure(self):
        self.filename = self.GetParameter("Filename")
        self.frame_number = 0

    def Physics(self, frame):

        logging.log_info('Storing MC truth information to a text file.')

        # Get the primary
        mctree = frame['I3MCTree']
        true_particle = mctree.primaries[0]
        
        # Get the daughters (currently only hadronic blobs)
        daughters = mctree.get_daughters(true_particle)

        assert(len(daughters) == 2)

        for daughter in daughters:
            # print(daughter)
            if daughter.type == dataclasses.I3Particle.Hadrons:
                # print("storing first cascade")
                casc_0_true = daughter
            else:
                # print("storing second cascade")
                # casc_1_true = daughter
                # MODIFICATION (for 190605 set, which has two hadron objects and intermediate HNL particle)
                casc_1_true = mctree.get_daughters(daughter)[0]
                
        # open true mc file
        # true_mc_file = open("{:s}_{:03d}_true_mc_info.txt".format(self.filename, self.frame_number), "w+")
        true_mc_file = open("{:s}_true_mc_info.txt".format(self.filename), "w+")

        true_mc_file.write("MC primary \n")
        true_mc_file.write(str(true_particle))
        true_mc_file.write("\n MC cascade 0 \n")
        true_mc_file.write(str(casc_0_true))
        true_mc_file.write("\n MC cascade 1 \n")
        true_mc_file.write(str(casc_1_true))
        true_mc_file.write("\n Scan Center \n")
        true_mc_file.write(str(frame[Seed]))
        true_mc_file.write("\n MC Truth \n")
        true_mc_file.write(str(frame[MCTruth]))

        true_mc_file.close()

        self.PushFrame(frame)
        self.frame_number += 1


def calculate_decayL(frame):
    mctree = frame['I3MCTree']
    true_particle = mctree.primaries[0]
    daughters = mctree.get_daughters(true_particle)
    # There are only ever 2 cascades in this simulation
    # First index should correspond to the first cascade at HNL production
    # and second should be the HNL decay vertex.      
    
    assert(len(daughters) == 2)
    
    for daughter in daughters:
        if daughter.type == dataclasses.I3Particle.Hadrons:
            casc_0_true = daughter
        else:
            # casc_1_true = daughter
            # MODIFICATION (for 190605 set, which has two hadron objects and intermediate HNL particle)
            casc_1_true = mctree.get_daughters(daughter)[0]
            
    decay_length = phys_services.I3Calculator.distance(casc_0_true,casc_1_true)/I3Units.m
    # Store in the frame for easy resolution checks later
    frame.Put("DecayLength",dataclasses.I3Double(decay_length))
    return True


def CreateTruthSeed(frame, seedname='TruthSeedParticle'):
    
    '''
    Create seed particle from primary truth information
    '''
    logging.log_info('Creating truth seed.')
    
    mctree = frame['I3MCTree']
    true_particle = mctree.primaries[0]

    # set truth values
    seed_particle = true_particle
    # set true decay length
    seed_particle.length = frame['DecayLength'].value
    seed_particle.fit_status = I3Particle.OK

    frame.Put(seedname, seed_particle)
    return True


##########      Main script     ##########

# GCDFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/SterileNeutrino/IC86/HighEnergy/HNL/MC/GCD/GeoCalibDetectorStatus_IC86.AVG_Pass2_SF0.99.i3.gz"
GCDFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz"

# M_1GeV is now 190600
# INFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/SterileNeutrino/IC86/HighEnergy/HNL/MC/M_1GeV/Ares/IC86.AVG/L7/domeff_0.97/00001-01000/L7_00_11_00056.i3.zst"
# INFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/double_bang_signal/00056_llh_scans/L7_00_11_00056_full_timed_fit_truth_seed.i3.zst"

# 190605 (using split fit v2 output files)
# INFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/190605/split_fit_v2/files/L7_00_11_01191.i3.zst"
# INFILE = "/lustre/fs22/group/icecube/lfischer/data/HNL/190605/split_fit_v2/files/L7_00_11_01791.i3.zst"

# INFILE = "/afs/ifh.de/user/l/lfischer/scratch/data/190605/taupede/files/L7_00_11_01191.i3.zst"
# INFILE = "/afs/ifh.de/user/l/lfischer/scratch/data/190605/taupede/files/L7_00_11_01791.i3.zst"
INFILE = "/afs/ifh.de/user/l/lfischer/scratch/data/190605/taupede/files/L7_00_11_02153.i3.zst"


# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2020/10_taupede_llh_scans/gulliview/selected_around_best_fit/L7_00_11_00056_around_truth_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2020/10_taupede_llh_scans/gulliview/selected_around_best_fit/L7_00_11_00056_around_best_fit'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/truth_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/truth_seed_inice_pulses'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2/truth_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2/truth_seed_split_inice_pulses'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/retro_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/retro_seed_inice_pulses'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2/retro_seed_split_inice_pulses'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/taupeded_reco_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1/taupeded_reco_seed_inice_pulses'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2/taupede_reco_seed_split_inice_pulses'


# # FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1_new/around_retro_seed'
# # FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1_new/around_mc_truth'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_1_new/around_taupede_bf'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2_new/around_retro_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2_new/around_mc_truth'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_2_new/around_taupede_bf'

# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_4/around_retro_seed'
# FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_4/around_mc_truth'
FILENAME_BASE = '/afs/ifh.de/user/l/lfischer/scratch/analysis_phd/plots_all/2021/01_taupede_signal/gulliview/event_4/around_taupede_bf'


# Pulses='SRTTWOfflinePulsesDC'
Pulses='SplitInIcePulses'

# Choose Scan Center
# Seed = 'RetroSeedParticle'
# Seed = 'TruthSeedParticle'
Seed = 'TaupedeFit'

# Set MC Truth Particle
MCTruth = 'TruthSeedParticle'

BadDOMs=["BadDomsList"]

# first set of steps
StepT=15
StepD=6
StepZenith=10
StepAzimuth=10
StepL=15

# # fine steps
# StepT=5
# StepD=1
# StepZenith=0.5
# StepAzimuth=0.5
# StepL=1

# # better steps
# StepT=5
# StepD=2
# StepZenith=3
# StepAzimuth=3
# StepL=6

inputfiles = [GCDFILE, INFILE]

tray = I3Tray()
tray.context["I3FileStager"] = icecube.dataio.get_stagers()

tray.Add("I3Reader", 'reader', filenamelist=inputfiles)

# Needed if only single test event should be selected
tray.AddModule(SelectEvent, "SelectEvent", Streams=[icetray.I3Frame.Physics])

# # Run the LLH-Scans on all contained events
# tray.AddModule(ContainmentFlags, "ContainmentFlags", Streams=[icetray.I3Frame.Physics])

tray.AddModule(calculate_decayL, "calculate_decayL", Streams=[icetray.I3Frame.Physics])

tray.AddModule(CreateTruthSeed, 'CreateTruthSeed', Streams=[icetray.I3Frame.Physics])

tray.Add(StoreMCTruth,
         Filename=FILENAME_BASE
         )

# old tables
# table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/ems_spice1_z20_a10.%s.fits')
# cascade_service = photonics_service.I3PhotoSplineService(table_base % 'abs', table_base % 'prob', 0, maxRadius=480.)

# new tables
table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
                                                         timingtable = table_base.format('prob'),
                                                         effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
                                                         timingSigma = 0,
                                                         maxRadius=400.
                                                        )

# add par
tray.AddService('TauMillipedeParametrizationFactory', 'paraname',
StepX=StepD*I3Units.m,
StepY=StepD*I3Units.m,
StepZ=StepD*I3Units.m,
StepT=StepT*I3Units.ns,
StepLinL=StepL*I3Units.m,
StepZenith=StepZenith*I3Units.degree,
StepAzimuth=StepAzimuth*I3Units.degree,
)

# add seeder
seed_kwargs = dict(InputReadout=Pulses, TimeShiftType="TNone", PositionShiftType="None")
seed_kwargs['FirstGuess'] = Seed
tray.AddService("I3BasicSeedServiceFactory", 'seeder', **seed_kwargs)

# add mc truth
seed_kwargs = dict(InputReadout=Pulses, TimeShiftType="TNone", PositionShiftType="None")
seed_kwargs['FirstGuess'] = MCTruth
tray.AddService("I3BasicSeedServiceFactory", 'mctruth', **seed_kwargs)

# add llh
millipede_config = dict(CascadePhotonicsService=cascade_service, Pulses=Pulses, ExcludedDOMs=list(set(['CalibrationErrata', 'SaturationWindows'] + BadDOMs)))
tray.AddService('MillipedeLikelihoodFactory', 'llhname', **millipede_config)

tray.Add(icecube.gulliver_modules.gulliview.GulliView,
            SeedService='seeder',
            MCTruth='mctruth',
            Parametrization='paraname',
            LogLikelihood='llhname',
            Filename=FILENAME_BASE)

tray.Execute()
tray.Finish()

t1 = time.time()
logging.log_info('Total time it took: {:.3f} s'.format(t1-t0))


##########      The End         ##########
