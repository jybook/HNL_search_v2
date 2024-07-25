#!/usr/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
#METAPROJECT /data/user/lfischer/software/icetray_main_real/build/
# this is a fresh icetray metaproject


# segmented-spline-reco metaproject with photospline release v2.1.0
###  #!/usr/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
###  #METAPROJECT /data/user/lfischer/software/icetray_main/build/


# tianlus metaproject
###  #!/usr/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/icetray-start
###  #METAPROJECT /data/user/tyuan/temp/tarballs/icetray.main.r323f33be.Linux-x86_64.gcc-9.3.0


# oscnext hnl metaproject
###  #!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
###  #METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/


### The Taupede (+MuMillipede) reconstruction script ###

#
# Author: Leander Fischer
#

import time
t0 = time.time()

from optparse import OptionParser

import os
import numpy as np
from icecube import icetray, dataclasses, dataio, photonics_service, phys_services, recclasses
from I3Tray import I3Tray
from icecube.dataclasses import I3Particle
from icecube.icetray import I3Bool

from icecube.millipede import TaupedeFit

import icecube.icetray.i3logging as logging

from copy import copy
from icecube.hdfwriter import I3HDFWriter


########## Functions needed to store as hdf5 ##########

HDF5_KEYS = [
    # EventProperties
    "decay_channel",
    "distance",
    "distanceMax",
    "distanceMin",
    "finalStateX",
    "finalStateY",
    "final_state_particle0",
    "final_state_particle1",
    "primary_type",
    "lifetime",
    "mHNL",
    "outgoing_neutrino_energy",
    "totalEnergy",
    "physical",
    "total_column_depth",
    # I3MCTree
    # true HNL variables
    "HNL_true_x",
    "HNL_true_y",
    "HNL_true_z",
    "HNL_true_energy",
    "HNL_true_zenith",
    "HNL_true_azimuth",
    "HNL_true_time",
    # true primary variables
    "true_x",
    "true_y",
    "true_z",
    "true_energy",
    "true_zenith",
    "true_azimuth",
    "true_time",
    # true first (DIS) cascade variables
    "casc0_true_x",
    "casc0_true_y",
    "casc0_true_z",
    "casc0_true_energy",
    "casc0_true_zenith",
    "casc0_true_azimuth",
    "casc0_true_time",
    # true second (HNL decay) cascade variables
    "casc1_true_x",
    "casc1_true_y",
    "casc1_true_z",
    "casc1_true_energy",
    "casc1_true_zenith",
    "casc1_true_azimuth",
    "casc1_true_time",
    "nan_decay_energy",
    # weights
    "LeptonInjectorWeight",
    "LifetimeWeight_1e-03",
    "OneWeight",
    "ReferenceWeight_1e-03",
    # total true quantities
    'true_decay_length',
    'true_total_energy',
    # reco first cascade variables
    "casc0_reco_x",
    "casc0_reco_y",
    "casc0_reco_z",
    "casc0_reco_energy",
    "casc0_reco_zenith",
    "casc0_reco_azimuth",
    "casc0_reco_time",
    # reco second cascade variables
    "casc1_reco_x",
    "casc1_reco_y",
    "casc1_reco_z",
    "casc1_reco_energy",
    "casc1_reco_zenith",
    "casc1_reco_azimuth",
    "casc1_reco_time",
    # total reco quantities
    'reco_decay_length',
    'reco_total_energy',
    # taupede fit parameters
    'taupede_logl',
    'taupede_rlogl',
    'taupede_ndof',
    'taupede_qtotal',
    'taupede_predicted_qtotal',
    'taupede_squared_residuals',
    'taupede_chi_squared',
    'taupede_chi_squared_dof',
    'taupede_nmini',
    # retro seed
    'retro_seed_x',
    'retro_seed_y',
    'retro_seed_z',
    'retro_seed_zenith',
    'retro_seed_azimuth',
    'retro_seed_length',
    'retro_seed_energy',
    'retro_seed_time',
    # flercnn seed
    'flercnn_seed_x',
    'flercnn_seed_y',
    'flercnn_seed_z',
    'flercnn_seed_zenith',
    'flercnn_seed_azimuth',
    'flercnn_seed_length',
    'flercnn_seed_energy',
    'flercnn_seed_time',
    # additional length seed options
    'L5_DirectHitsC_DirTrackLength',
    'FiniteRecoFit_length',
    # PID usable variables
    'L7_PIDClassifier_FullSky_ProbTrack',
    'L7_PIDClassifier_Upgoing_ProbTrack',
    'retro_crs_prefit__zero_dllh_mean',
    'retro_crs_prefit__max_llh',
    # MuMillipede variables
    'E_tot',
    'E_cascade0_5m',
    'E_cascade1_5m',
    'E_cascade0_10m',
    'E_cascade1_10m',
    'E_cascade0_20m',
    'E_cascade1_20m',
    'E_cascade0_40m',
    'E_cascade1_40m',
]

# remove duplicates from the HDF5 keys
HDF5_KEYS = list(
    set(HDF5_KEYS)  # note that this does NOT preserve order, but doesn't matter
)

# true variables - model independent simulation case
def write_true_data_to_keys(frame):
    """
    Function to write true cascade information, energy and decay length to frame
    """
    cascades = frame["I3MCTree"].get_daughters(frame["I3MCTree"].primaries[0])

    # add true cascades, energy and decay length to the frame
    frame["true_decay_length"] = dataclasses.I3Double(
        phys_services.I3Calculator.distance(cascades[0], cascades[1])
    )
    frame["true_total_energy"] = dataclasses.I3Double(
        cascades[0].energy + cascades[1].energy
    )
    for count, cascade in enumerate(cascades):
        frame["casc{}_true_x".format(count)] = dataclasses.I3Double(cascade.pos.x)
        frame["casc{}_true_y".format(count)] = dataclasses.I3Double(cascade.pos.y)
        frame["casc{}_true_z".format(count)] = dataclasses.I3Double(cascade.pos.z)
        frame["casc{}_true_energy".format(count)] = dataclasses.I3Double(cascade.energy)
        frame["casc{}_true_zenith".format(count)] = dataclasses.I3Double(cascade.dir.zenith)
        frame["casc{}_true_azimuth".format(count)] = dataclasses.I3Double(cascade.dir.azimuth)
        frame["casc{}_true_time".format(count)] = dataclasses.I3Double(cascade.time)

    return True

# true variables - model specific HNL simulation case
def store_mc_true_variables(frame):
    if not frame.Has("I3MCTree"):
        return False
    mctree = frame["I3MCTree"]
    p_true = mctree.primaries[0]

    p_daughters = mctree.get_daughters(p_true)
    assert len(p_daughters) == 2

    for p_daughter in p_daughters:
        if p_daughter.type == dataclasses.I3Particle.Hadrons:
            casc_0_true = p_daughter
        else:
            hnl_true = p_daughter

    hnl_daughters = mctree.get_daughters(hnl_true)
    assert len(hnl_daughters) > 0

    nan_energy = False

    for count_hnl_daughters, hnl_daughter in enumerate(hnl_daughters):
        if not count_hnl_daughters:
            casc_1_true = hnl_daughter
        else:
            assert casc_1_true.pos == hnl_daughter.pos
            if np.isnan(hnl_daughter.energy):
                nan_energy = True
            casc_1_true.energy = casc_1_true.energy + hnl_daughter.energy

    frame["true_x"] = dataclasses.I3Double(p_true.pos.x)
    frame["true_y"] = dataclasses.I3Double(p_true.pos.y)
    frame["true_z"] = dataclasses.I3Double(p_true.pos.z)
    frame["true_energy"] = dataclasses.I3Double(p_true.energy)
    frame["true_zenith"] = dataclasses.I3Double(p_true.dir.zenith)
    frame["true_azimuth"] = dataclasses.I3Double(p_true.dir.azimuth)
    frame["true_time"] = dataclasses.I3Double(p_true.time)

    frame["HNL_true_x"] = dataclasses.I3Double(hnl_true.pos.x)
    frame["HNL_true_y"] = dataclasses.I3Double(hnl_true.pos.y)
    frame["HNL_true_z"] = dataclasses.I3Double(hnl_true.pos.z)
    frame["HNL_true_energy"] = dataclasses.I3Double(hnl_true.energy)
    frame["HNL_true_zenith"] = dataclasses.I3Double(hnl_true.dir.zenith)
    frame["HNL_true_azimuth"] = dataclasses.I3Double(hnl_true.dir.azimuth)
    frame["HNL_true_time"] = dataclasses.I3Double(hnl_true.time)

    frame["casc0_true_x"] = dataclasses.I3Double(casc_0_true.pos.x)
    frame["casc0_true_y"] = dataclasses.I3Double(casc_0_true.pos.y)
    frame["casc0_true_z"] = dataclasses.I3Double(casc_0_true.pos.z)
    frame["casc0_true_energy"] = dataclasses.I3Double(casc_0_true.energy)
    frame["casc0_true_zenith"] = dataclasses.I3Double(casc_0_true.dir.zenith)
    frame["casc0_true_azimuth"] = dataclasses.I3Double(casc_0_true.dir.azimuth)
    frame["casc0_true_time"] = dataclasses.I3Double(casc_0_true.time)

    frame["casc1_true_x"] = dataclasses.I3Double(casc_1_true.pos.x)
    frame["casc1_true_y"] = dataclasses.I3Double(casc_1_true.pos.y)
    frame["casc1_true_z"] = dataclasses.I3Double(casc_1_true.pos.z)
    frame["casc1_true_energy"] = dataclasses.I3Double(casc_1_true.energy)
    frame["casc1_true_zenith"] = dataclasses.I3Double(casc_1_true.dir.zenith)
    frame["casc1_true_azimuth"] = dataclasses.I3Double(casc_1_true.dir.azimuth)
    frame["casc1_true_time"] = dataclasses.I3Double(casc_1_true.time)

    frame["nan_decay_energy"] = icetray.I3Bool(nan_energy)

    return True

# function to get EventProperties (LI) - this won't run with the standard icetray software
def store_LI_event_properties(frame):
    if not frame.Has("EventProperties"):
        return False
    event_properties = frame["EventProperties"]

    frame["mHNL"] = dataclasses.I3Double(event_properties.mHNL)
    frame["distance"] = dataclasses.I3Double(event_properties.distance)
    frame["distanceMin"] = dataclasses.I3Double(event_properties.distanceMin)
    frame["distanceMax"] = dataclasses.I3Double(event_properties.distanceMax)
    frame["lifetime"] = dataclasses.I3Double(event_properties.lifetime)
    frame["totalEnergy"] = dataclasses.I3Double(event_properties.totalEnergy)
    frame["finalStateX"] = dataclasses.I3Double(event_properties.finalStateX)
    frame["finalStateY"] = dataclasses.I3Double(event_properties.finalStateY)
    frame['primary_type'] =  dataclasses.I3Double(event_properties.initialType)
    frame['final_state_particle0'] =  dataclasses.I3Double(event_properties.finalType1)
    frame['final_state_particle1'] =  dataclasses.I3Double(event_properties.finalType2)
    frame['total_column_depth'] =  dataclasses.I3Double(event_properties.totalColumnDepth)
    frame["decay_channel"] = dataclasses.I3Double(event_properties.decay_channel)
    frame["outgoing_neutrino_energy"] = dataclasses.I3Double(event_properties.outgoingNeutrinoEnergy)
    frame["physical"] = dataclasses.I3Double(event_properties.physical)

    return True

# function to get weights
def store_weights(frame):
    if not frame.Has("I3MCWeightDict"):
        return False
    weight_dict = frame["I3MCWeightDict"]

    frame["LeptonInjectorWeight"] = dataclasses.I3Double(weight_dict['LeptonInjectorWeight'])
    frame["LifetimeWeight_1e-03"] = dataclasses.I3Double(weight_dict['LifetimeWeight_1e-03'])
    frame["OneWeight"] = dataclasses.I3Double(weight_dict['OneWeight'])
    frame["ReferenceWeight_1e-03"] = dataclasses.I3Double(weight_dict['ReferenceWeight_1e-03'])

    return True

def write_reco_data_to_keys(frame, fitname='TaupedeFit'):
    """
    Function to write reco cascade information, energy and decay length to frame
    """
    fit_cascades_key = fitname + 'Particles'

    if frame.Has(fit_cascades_key):
        cascades = frame[fit_cascades_key]
               
        # add reco cascades, energy and decay length to the frame
        frame["reco_decay_length"] = dataclasses.I3Double(
            phys_services.I3Calculator.distance(cascades[0], cascades[1])
        )
        frame["reco_total_energy"] = dataclasses.I3Double(
            cascades[0].energy + cascades[1].energy
        )
        for count, cascade in enumerate(cascades):
            frame["casc{}_reco_x".format(count)] = dataclasses.I3Double(cascade.pos.x)
            frame["casc{}_reco_y".format(count)] = dataclasses.I3Double(cascade.pos.y)
            frame["casc{}_reco_z".format(count)] = dataclasses.I3Double(cascade.pos.z)
            frame["casc{}_reco_energy".format(count)] = dataclasses.I3Double(cascade.energy)
            frame["casc{}_reco_zenith".format(count)] = dataclasses.I3Double(cascade.dir.zenith)
            frame["casc{}_reco_azimuth".format(count)] = dataclasses.I3Double(cascade.dir.azimuth)
            frame["casc{}_reco_time".format(count)] = dataclasses.I3Double(cascade.time)

    # add some more reco data (from previous fits)
    if frame.Has('L5_DirectHitsC'):
        frame["L5_DirectHitsC_DirTrackLength"] = dataclasses.I3Double(frame['L5_DirectHitsC'].dir_track_length)
    if frame.Has('FiniteRecoFit'):
        frame["FiniteRecoFit_length"] = dataclasses.I3Double(frame['FiniteRecoFit'].length)
    if frame.Has('retro_crs_prefit__zero_dllh'):
        frame['retro_crs_prefit__zero_dllh_mean'] = dataclasses.I3Double(frame['retro_crs_prefit__zero_dllh']['mean'])

    return True

def write_taupede_fitparameters(frame, fitname='TaupedeFit'):
    """
    Function to write taupede goodness of fit parameters
    """

    fit_params_key = fitname + 'FitParams'
    if frame.Has(fit_params_key):
        fit_params = frame[fit_params_key]
        frame['taupede_logl'] = dataclasses.I3Double(fit_params.logl)
        frame['taupede_rlogl'] = dataclasses.I3Double(fit_params.rlogl)
        frame['taupede_ndof'] = dataclasses.I3Double(fit_params.ndof)
        frame['taupede_qtotal'] = dataclasses.I3Double(fit_params.qtotal)
        frame['taupede_predicted_qtotal'] = dataclasses.I3Double(fit_params.predicted_qtotal)
        frame['taupede_squared_residuals'] = dataclasses.I3Double(fit_params.squared_residuals)
        frame['taupede_chi_squared'] = dataclasses.I3Double(fit_params.chi_squared)
        frame['taupede_chi_squared_dof'] = dataclasses.I3Double(fit_params.chi_squared_dof)
        frame['taupede_nmini'] = dataclasses.I3Double(fit_params.nmini)

    return True

def extract_MuMillipede_energies(frame, MuMillipede_fitname='MuMillipedeFit', Taupede_fitname='TaupedeFit'):
    """
    Function to extract energies around the fitted cascades from the infinite track fit
    """

    fit_cascades_key = Taupede_fitname + 'Particles'
    if frame.Has(fit_cascades_key):
        if frame.Has(MuMillipede_fitname):
            trackvector = frame[MuMillipede_fitname]

            E_tot = 0.0
            E_cascade0_5m = 0.0
            E_cascade1_5m = 0.0
            E_cascade0_10m = 0.0
            E_cascade1_10m = 0.0
            E_cascade0_20m = 0.0
            E_cascade1_20m = 0.0
            E_cascade0_40m = 0.0
            E_cascade1_40m = 0.0

            for sec in trackvector:
                E_tot += sec.energy
                d0 = np.sqrt(
                    np.power(sec.pos.x - frame['casc0_reco_x'].value, 2)
                    + np.power(sec.pos.y - frame['casc0_reco_y'].value, 2)
                    + np.power(sec.pos.z - frame['casc0_reco_z'].value, 2)
                )
                d1 = np.sqrt(
                    np.power(sec.pos.x - frame['casc1_reco_x'].value, 2)
                    + np.power(sec.pos.y - frame['casc1_reco_y'].value, 2)
                    + np.power(sec.pos.z - frame['casc1_reco_z'].value, 2)
                )
                # 5 m
                if d0 < 5 and d1 > 5:
                    E_cascade0_5m += sec.energy
                elif d0 > 5 and d1 < 5:
                    E_cascade1_5m += sec.energy
                elif d0 < 5 and d1 < 5:
                    if d0 < d1:
                        E_cascade0_5m += sec.energy
                    else:
                        E_cascade1_5m += sec.energy
                # 10 m
                if d0 < 10 and d1 > 10:
                    E_cascade0_10m += sec.energy
                elif d0 > 10 and d1 < 10:
                    E_cascade1_10m += sec.energy
                elif d0 < 10 and d1 < 10:
                    if d0 < d1:
                        E_cascade0_10m += sec.energy
                    else:
                        E_cascade1_10m += sec.energy
                # 20 m
                if d0 < 20 and d1 > 20:
                    E_cascade0_20m += sec.energy
                elif d0 > 20 and d1 < 20:
                    E_cascade1_20m += sec.energy
                elif d0 < 20 and d1 < 20:
                    if d0 < d1:
                        E_cascade0_20m += sec.energy
                    else:
                        E_cascade1_20m += sec.energy
                # 40 m
                if d0 < 40 and d1 > 40:
                    E_cascade0_40m += sec.energy
                elif d0 > 40 and d1 < 40:
                    E_cascade1_40m += sec.energy
                elif d0 < 40 and d1 < 40:
                    if d0 < d1:
                        E_cascade0_40m += sec.energy
                    else:
                        E_cascade1_40m += sec.energy

            frame['E_tot'] = dataclasses.I3Double(E_tot)
            frame['E_cascade0_5m'] = dataclasses.I3Double(E_cascade0_5m)
            frame['E_cascade1_5m'] = dataclasses.I3Double(E_cascade1_5m)
            frame['E_cascade0_10m'] = dataclasses.I3Double(E_cascade0_10m)
            frame['E_cascade1_10m'] = dataclasses.I3Double(E_cascade1_10m)
            frame['E_cascade0_20m'] = dataclasses.I3Double(E_cascade0_20m)
            frame['E_cascade1_20m'] = dataclasses.I3Double(E_cascade1_20m)
            frame['E_cascade0_40m'] = dataclasses.I3Double(E_cascade0_40m)
            frame['E_cascade1_40m'] = dataclasses.I3Double(E_cascade1_40m)

            return True
        return False
    return False


########## End


### Start - Re-define HighEnergyExclusions ###
 
@icetray.traysegment
def HighEnergyExclusions(
    tray,
    name,
    Pulses,
    ExcludeDeepCore=None,
    ExcludeSaturatedDOMs=None,
    ExcludeBrightDOMs='BrightDOMs',
    BrightDOMThreshold=30,
    SaturationWindows='SaturationWindows',
    BadDomsList='BadDomsList',
    CalibrationErrata='CalibrationErrata'
    ):
    """
    Work around systematic errors in the modelling of the detector response by removing certain classes
    of DOMs from consideration that would otherwise over-contribute to Millipede likelihoods for events
    above a few hundred TeV.

    The options beginning with "Exclude" may be set to None or False to disable the relevant exclusion.

    :param Pulses: the name of the pulse map to be used for reconstruction
    :param ExcludeDeepCore: remove DeepCore strings from consideration
    :param ExcludeSaturatedDOMs: exclude saturated DOMs entirely, not just during the times when their output current is above the linearity limit
    :param ExcludeBrightDOMs: exclude DOMs that collect a total charge a factor greater than the mean charge
    :param BrightDOMThreshold: threshold factor for bright DOMs
    :param BadDomsList: list of DOMs that can't produce useful data
    :param SaturationWindows: times during which PMTs were nonlinear

    :returns: a list of exclusions that can be passed to Millipede modules
    """

    def log(message):
        icetray.logging.log_info(message, unit="MillipedeHighEnergyExclusions")

    exclusions = [CalibrationErrata, BadDomsList]
    if ExcludeDeepCore:
        log("Excluding DeepCore DOMs, are you sure?")
        def DeepCoreFlagger(frame):
            dc_strings = set(range(79,87))
            dc_doms = dataclasses.I3VectorOMKey()
            for om in frame['I3Geometry'].omgeo.keys():
                if om.string in dc_strings:
                    dc_doms.append(om)
            frame[ExcludeDeepCore] = dc_doms
        tray.AddModule(DeepCoreFlagger, name+'DeepCoreFlagger', Streams=[icetray.I3Frame.Geometry])
        exclusions.append(ExcludeDeepCore)
	
    if ExcludeSaturatedDOMs:
        def SaturationExtender(frame):
            if SaturationWindows in frame:
                frame[ExcludeSaturatedDOMs] = dataclasses.I3VectorOMKey(frame[SaturationWindows].keys())
                log("Excluding %d saturated DOMs: %s" % (len(frame[ExcludeSaturatedDOMs]), list(frame[ExcludeSaturatedDOMs])))
        tray.AddModule(SaturationExtender, name+'SaturationExtender', Streams=[icetray.I3Frame.Physics])
        exclusions.append(ExcludeSaturatedDOMs)
    else:
        exclusions.append(SaturationWindows)
	
    if ExcludeBrightDOMs:
        def BrightDOMFlagger(frame):
            bright_doms = dataclasses.I3VectorOMKey()
            pmap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, Pulses)

            # now using median, wow so sophisticated
            q_median = np.median( [ sum(p.charge for p in pulses) for pulses in pmap.values() ] )
            max_q = BrightDOMThreshold*q_median

            for om, pulses in pmap.items():
                q = sum(p.charge for p in pulses)
                if q > max_q:
                    bright_doms.append(om)
            if len(bright_doms):
                frame[ExcludeBrightDOMs] = bright_doms
                log("Excluding %d bright DOMs: %s" % (len(bright_doms), list(bright_doms)))
                
        tray.AddModule(BrightDOMFlagger, name+'BrightDOMFlagger')
        exclusions.append(ExcludeBrightDOMs)

    return exclusions

### End - Re-define HighEnergyExclusions ###


### Start - Flip reco cascade functions ###

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
    Function to correct for negative reconstructed decay lengths/flipped cascade order.
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

### End - Flip reco cascade functions ###


### Start - Create seeds functions ###

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

    frame["retro_seed_x"] = dataclasses.I3Double(seed_particle.pos.x)
    frame["retro_seed_y"] = dataclasses.I3Double(seed_particle.pos.y)
    frame["retro_seed_z"] = dataclasses.I3Double(seed_particle.pos.z)
    
    # set direction
    seed_particle.dir = dataclasses.I3Direction(
        frame[retro_key_base+'zenith'].value,
        frame[retro_key_base+'azimuth'].value,
    )

    frame["retro_seed_zenith"] = dataclasses.I3Double(seed_particle.dir.zenith)
    frame["retro_seed_azimuth"] = dataclasses.I3Double(seed_particle.dir.azimuth)

    # set track length
    seed_particle.length = frame[retro_key_base+'track_length'].value

    frame["retro_seed_length"] = dataclasses.I3Double(seed_particle.length)
    
    # set energy as em cascade energy
    seed_particle.energy = frame[retro_key_base+'em_cascade_energy'].value

    frame["retro_seed_energy"] = dataclasses.I3Double(seed_particle.energy)

    # set time
    seed_particle.time = frame[retro_key_base+'time'].value

    frame["retro_seed_time"] = dataclasses.I3Double(seed_particle.time)

    # set fit status
    seed_particle.shape = I3Particle.Primary
    seed_particle.fit_status = I3Particle.OK

    # set particle type
    seed_particle.type = I3Particle.NuTau
    seed_particle.pdg_encoding = 16

    seed_particle.speed = dataclasses.I3Constants.c

    frame.Put(seedname, seed_particle)
    return True

def CreateFlercnnSeed(frame, seedname='FlercnnSeedParticle'):
    '''
    Create seed particle from flercnn best fit particle with em cascade energy
    '''
    logging.log_info('Creating flercnn reco seed.')

    flercnn_key_base = 'FLERCNN_'

    seed_particle = I3Particle()
    # set vertex
    seed_particle.pos = dataclasses.I3Position(
        frame[flercnn_key_base+'vertex_x'].value,
        frame[flercnn_key_base+'vertex_y'].value,
        frame[flercnn_key_base+'vertex_z'].value,
    )

    frame["flercnn_seed_x"] = dataclasses.I3Double(seed_particle.pos.x)
    frame["flercnn_seed_y"] = dataclasses.I3Double(seed_particle.pos.y)
    frame["flercnn_seed_z"] = dataclasses.I3Double(seed_particle.pos.z)
    
    tmp_azimuth = 0.0  # set to zero and overwrite with FiniteRecoFit
    tmp_time = 10000.0  # set default value
    tmp_length = 100.0  # set default value

    if frame.Has('FiniteRecoFit'):
        FiniteRecoFit_particle = frame['FiniteRecoFit']
        if FiniteRecoFit_particle.fit_status == I3Particle.OK:
            logging.log_info('FiniteRecoFit is present and valid, use it to seed azimuth and time.')
            # set azimuth
            tmp_azimuth = FiniteRecoFit_particle.dir.azimuth
            # set time
            tmp_time = FiniteRecoFit_particle.time
            # set track length
            tmp_length = FiniteRecoFit_particle.length
        else:
            logging.log_info('FiniteRecoFit is present but NOT valid, use default azimuth/time/length.')

    if frame.Has('OnlineL2_SplineMPE_DirectHitsC'):
        logging.log_info('OnlineL2_SplineMPE_DirectHitsC is present, use it to seed track length.')
        # overwrite track length
        tmp_length = frame['OnlineL2_SplineMPE_DirectHitsC'].dir_track_length

    # set the values for the seed particle from FiniteRecoFit/SplineMPE/default
    # set time
    seed_particle.time = tmp_time
    # set track length
    seed_particle.length = tmp_length
    
    # set zenith (and dummy azimuth)
    seed_particle.dir = dataclasses.I3Direction(
        frame[flercnn_key_base+'zenith'].value,
        tmp_azimuth,  # set to zero and overwrite with FiniteRecoFit
    )

    frame["flercnn_seed_zenith"] = dataclasses.I3Double(seed_particle.dir.zenith)
    frame["flercnn_seed_azimuth"] = dataclasses.I3Double(seed_particle.dir.azimuth)
    frame["flercnn_seed_time"] = dataclasses.I3Double(seed_particle.time)
    frame["flercnn_seed_length"] = dataclasses.I3Double(seed_particle.length)
    
    # set energy as em cascade energy
    seed_particle.energy = frame[flercnn_key_base+'energy'].value

    frame["flercnn_seed_energy"] = dataclasses.I3Double(seed_particle.energy)

    # set fit status
    seed_particle.shape = I3Particle.Primary
    seed_particle.fit_status = I3Particle.OK

    # set particle type
    seed_particle.type = I3Particle.NuTau
    seed_particle.pdg_encoding = 16

    seed_particle.speed = dataclasses.I3Constants.c

    frame.Put(seedname, seed_particle)
    return True

def CreateTruthSeed(frame, randomservice, seedname='TruthSeedParticle'):
    '''
    Create seed particle from MC truth
    '''
    logging.log_info('Creating MC truth seed.')

    cascade_0 = frame["I3MCTree"].get_daughters(frame["I3MCTree"].primaries[0])[0]
    true_decay_length = frame['true_decay_length'].value
    total_energy = frame['true_total_energy'].value

    seed_particle = I3Particle()
    seed_particle.pos = dataclasses.I3Position(
        randomservice.gaus(cascade_0.pos.x, 5.0),
        randomservice.gaus(cascade_0.pos.y, 5.0),
        randomservice.gaus(cascade_0.pos.z, 2.5),
    )
    seed_particle.dir = dataclasses.I3Direction(
        randomservice.gaus(cascade_0.dir.zenith, 0.1),
        randomservice.gaus(cascade_0.dir.azimuth, 0.2),
    )
    seed_particle.length = randomservice.gaus(true_decay_length, 5.0)
    seed_particle.energy = randomservice.gaus(total_energy, 0.3*total_energy)
    seed_particle.time = randomservice.gaus(cascade_0.time, 50.0)
    
    seed_particle.shape = I3Particle.Primary
    seed_particle.fit_status = I3Particle.OK
    seed_particle.type = I3Particle.NuTau
    seed_particle.pdg_encoding = 16
    seed_particle.speed = dataclasses.I3Constants.c

    frame.Put(seedname, seed_particle)
    return True

def CreateL2Seed(frame, newseedname):
    '''
    Create seed particle from FiniteRecoFit using OnlineL2_SplineMPE_DirectHitsC DirTrackLength as length
    '''
    logging.log_info('Creating L2 seed particle.')

    if frame.Has('FiniteRecoFit') and frame.Has('OnlineL2_SplineMPE_DirectHitsC'):
        logging.log_info('Both keys were there, using full L2 seed.')
        seed_particle = frame['FiniteRecoFit']
        seed_particle.length = frame['OnlineL2_SplineMPE_DirectHitsC'].dir_track_length
        seed_particle.shape = I3Particle.Primary
        seed_particle.fit_status = I3Particle.OK
        seed_particle.type = I3Particle.NuTau
        seed_particle.pdg_encoding = 16
        seed_particle.speed = dataclasses.I3Constants.c
        frame.Put('written_correct_L2_seed', I3Bool(True))
    else:
        logging.log_info('At least one of the input keys is missing, writing empty seed particle and returning False.')
        seed_particle = I3Particle()
        frame.Put('written_correct_L2_seed', I3Bool(False))
        return False

    frame.Put(newseedname, seed_particle)
    return True

# Add length to seed if it does not have one.
def AddLengthToSeed(frame, newseedname, seedparticle):
    '''
    Create seed particle from given name and add length if it does not have one.
    '''
    logging.log_info('Creating seed (making sure it has a length).')

    seed_particle = I3Particle(frame[seedparticle])

    if np.isnan(seed_particle.length):
        # set track length
        seed_particle.length = 100.

    frame.Put(newseedname, seed_particle)
    return True

### End - Create seeds functions ###


###########################
###   TAUPEDE
###########################

@icetray.traysegment
def TaupedeFitWrapper(tray, name, Seed, Iterations, PhotonsPerBin,
                      lengthbounds, **millipede_params):
    
    logging.log_info('Running Taupede fit.')

    """
    Chosen fit routine: Fit is split into length+energies fit and full fit,
    where the latter is using the result from the first as a seed.
    Length+energies fit are additionally seeded with 3 scaled values of the
    seed track length. Best performing is selected.
    """

    # 3-fold scaling
    scales = [0.5, 1., 1.5]
    # # more_but_coarser_length_seeds
    # # 5-fold scaling
    # scales = [0.2, 0.5, 1., 2., 5., 10.]

    def addbruteforceseeds(frame, Seed, Output):

        logging.log_info('Adding brute force seeds to the frame.')

        if frame.Has(Seed):
            taupedeseedbase = I3Particle(frame[Seed])
        else:
            logging.log_info("Using empty I3Particle because Seed is missing")
            taupedeseedbase = I3Particle()
        for scale in scales:
            taupedeseed = I3Particle(taupedeseedbase)
            taupedeseed.length = scale * taupedeseed.length
            frame.Put(Output+'_scale{:.01f}'.format(scale), taupedeseed)

    tray.Add(
        addbruteforceseeds,
        'AddBruteForceSeeds',
        Seed=Seed,
        Output='TaupedeBruteForceWrapperSeed',
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    for scale in scales:
        seedtag = 'TaupedeBruteForceWrapperSeed_scale{:.01f}'.format(scale)
        length_fit_tag = 'TaupedeBruteForceWrapperReco_scale{:.01f}'.format(scale)

        # length fit (not timed)
        tray.Add(
            TaupedeFit,
            length_fit_tag,
            Seed = seedtag,
            LengthBounds = lengthbounds,
            StepT=0,
            StepD=0,
            StepZenith=0,
            StepAzimuth=0,
            # initial step size for length
            StepL = 1,
            # # maybe an intermediate decay length is the way to go
            # StepL = 2,
            # # more_but_coarser_length_seeds
            # StepL = 5,
            Iterations = 5,
            PhotonsPerBin = -1,
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            **millipede_params
            )

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

        dtype = [('fittag', 'S60'), ('rlogl', float)]
        fitmap = np.array(rlogl_map, dtype=dtype)

        sortedfitmap = np.sort(fitmap, order='rlogl')
    
        if not len(sortedfitmap) == 0:
            bestfittag = sortedfitmap[0]['fittag']

            bestfitparticle = I3Particle(frame[bestfittag])
        else:
            bestfitparticle = I3Particle(frame[Seed])

        frame.Put(Output, bestfitparticle)

    tray.Add(
        findbestfit,
        'FindBestFit',
        Seed=Seed,
        Output='TaupedeBestLengthFit',
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )
    
    def CheckSeedIsValid(frame, SeedName='TaupedeBestLengthFit'):

        # logging.log_info('Checking whether the seed going into full Taupede fit is valid.')

        if frame.Has(SeedName):
            if frame[SeedName].fit_status == I3Particle.OK:
                pass
            else:
                logging.log_info("Seedparticle (length prefit) fit status is NOT OK, using original seed.")
                TaupedeSeed = I3Particle(frame[Seed])
                frame.Put(SeedName, TaupedeSeed)
        else:
            logging.log_info("Seedparticle (length prefit) is missing completely, using original seed.")
            TaupedeSeed = I3Particle(frame[Seed])
            frame.Put(SeedName, TaupedeSeed)

    tray.Add(
        CheckSeedIsValid,
        'CheckSeedIsValid',
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    recotag = name
  
    logging.log_info("Running full fit now:")

    # full fit
    tray.Add(
        TaupedeFit, 
        recotag,
        Seed = 'TaupedeBestLengthFit',
        LengthBounds = lengthbounds,
        StepL = 1,
        Iterations = Iterations,
        PhotonsPerBin = PhotonsPerBin,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        **millipede_params
        )

###########################

def main():

    """
    Perform reconstruction and pre-processing
    """

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--infile", dest="INFILE", help="Read input from INFILE (.i3{.gz/.zst} format) [default: %default]")
    parser.add_option("-o", "--outfile", dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz/.zst} format) [default: %default]")
    parser.add_option("-g", "--gcdfile",
                    default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz",
                    dest="GCDFILE", help="Read in GCD file [default: %default]")
    parser.add_option("-s", "--seedparticle", default='FlercnnSeedParticle',
                    dest="SEEDPARTICLE", help="Choose particle to use as seed for taupede, either real I3Particle or from ('RetroSeedParticle', 'FlercnnSeedParticle, 'TruthSeedParticle' or 'FiniteRecoFit+OnlineL2_SplineMPE_DirectHitsC_DirTrackLength') [default: %default]")
    parser.add_option("--randseed",dest="randseed",type="int",
                    help="Random seed(int) for the truth seed smearing. Will extract filenumber from filename if not set [default: %default]", default=None)
    parser.add_option("--BrightDOMThreshold",dest="BrightDOMThreshold",type="int",
                    help="DOMs with charge larger than this factor times the median will be excluded from fit [default: %default]", default=30)
    parser.add_option("-m","--icemodel", default="spice_lea", type='string',
                    dest="ICEMODEL",help="Ice model (spice_3.2.1, spice_lea etc.) [default: %default]")
    parser.add_option("-p", "--pulseseries", default='SplitInIcePulses',
                    dest="PULSESERIES", help="Choose pulse series to use [default: %default]")
    parser.add_option("-v", "--verbose", default=False, action="store_true",
                    dest="VERBOSE", help="Print info level output [default: %default]")
    parser.add_option("--hdf5", action="store_false",
                    dest="write_hdf5", default=True, help="Also store hdf5 file [default: %default]")
    parser.add_option("--real_sim", type="int",
                    dest="real_sim", default=1, help="Is this the real (model dependent) simulation? [default: %default]")
    (options,args) = parser.parse_args()
    if not options.INFILE:
        parser.error('Infile not specified.')
    if not options.GCDFILE:
        parser.error('GCDfile not specified.')

    write_hdf5  = options.write_hdf5
    real_sim    = bool(options.real_sim)
    randseed 	= options.randseed
    pulses      = options.PULSESERIES

    if options.VERBOSE:
        logging.set_level('INFO')
    else:
        logging.set_level('WARN')

    logging.log_info('Options: {}'.format(options))

    OUTPATH = os.path.split(options.OUTFILE)[0]
    if not os.path.exists(OUTPATH):
        logging.log_info('Creating output filepath: {}'.format(OUTPATH))
        os.makedirs(OUTPATH)

    if not randseed:
        randseed = int(options.INFILE.split('.i3')[0][-5:])
        logging.log_info('Using random seed: {}'.format(randseed))

    # set random seed (used for truth seed smearing)
    rndserv = phys_services.I3GSLRandomService(seed = randseed)


    ### START ###

    tray = I3Tray()

    tray.context['I3FileStager'] = dataio.get_stagers()

    tray.AddModule(
        'I3Reader',
        'reader',
        FilenameList=[options.GCDFILE, options.INFILE],
        )

    # # quick hack to only run on specific eventidID
    # def eventid_filter(frame):
    #     logging.log_info("Running eventid filter (eventid:{})".format(frame['I3EventHeader'].event_id))
    #     if not frame['I3EventHeader'].event_id == 1660:
    #         logging.log_info('Skipping eventid: {}'.format(frame['I3EventHeader'].event_id))
    #         return False

    # tray.AddModule(
    #     eventid_filter,
    #     'eventid_filter',
    #     If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit"
    #     )

    # get true simulation data (dependent on the simulation used)
    if real_sim:
        # add function to write true data as keys - model dependent case
        tray.AddModule(
            store_mc_true_variables,
            "WriteTrueDataToKeys_Model_Dependent",
            Streams=[icetray.I3Frame.Physics],
        )
    else:
        # add function to write true data as keys - model independent case
        tray.AddModule(
            write_true_data_to_keys,
            "WriteTrueDataToKeys_Model_Independent",
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )

    if options.SEEDPARTICLE == 'RetroSeedParticle':
        # set taupede seed particle from retro
        tray.AddModule(
            CreateRetroSeed,
            'CreateRetroSeed',
            Streams=[icetray.I3Frame.Physics],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        Seed = options.SEEDPARTICLE
    elif options.SEEDPARTICLE == 'FlercnnSeedParticle':
        # set taupede seed particle from flercnn
        tray.AddModule(
            CreateFlercnnSeed,
            'CreateFlercnnSeed',
            Streams=[icetray.I3Frame.Physics],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        Seed = options.SEEDPARTICLE
    elif options.SEEDPARTICLE == 'TruthSeedParticle':
        # set taupede seed particle from MC truth
        tray.AddModule(
            CreateTruthSeed,
            'CreateTruthSeed',
            randomservice = rndserv,
            Streams=[icetray.I3Frame.Physics],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        Seed = options.SEEDPARTICLE
    elif options.SEEDPARTICLE == 'FiniteRecoFit+OnlineL2_SplineMPE_DirectHitsC_DirTrackLength':
        # set taupede seed particle from FiniteRecoFit using OnlineL2_SplineMPE_DirectHitsC DirTrackLength as length
        newseedname='L2SeedParticle'
        tray.AddModule(
            CreateL2Seed,
            'CreateL2Seed',
            newseedname='L2SeedParticle',
            Streams=[icetray.I3Frame.Physics],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        Seed = newseedname
        HDF5_KEYS.append('written_correct_L2_seed')
    else:
        # add a decay length to the seed particle if it does not have one
        newseedname = 'TaupedeInitialSeed' + options.SEEDPARTICLE
        tray.AddModule(
            AddLengthToSeed,
            'AddLengthToSeed',
            newseedname=newseedname,
            seedname=options.SEEDPARTICLE,
            Streams=[icetray.I3Frame.Physics],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            )
        Seed = newseedname

    if options.ICEMODEL == 'spice_lea':
        logging.log_info('Using Spice Lea tables:')
        # add table base/cascade service (photonics tables including effective distande parametrization by Marcel Usner)
        # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
        table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
        cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
                                                                timingtable = table_base.format('prob'),
                                                                effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
                                                                timingSigma = 0,
                                                                maxRadius=400.,
                                                                )
    # if options.ICEMODEL == 'spice_lea_double_precision_1e-08':
    #     logging.log_info('Using Spice Lea tables:')
    #     # add table base/cascade service (photonics tables including effective distande parametrization by Marcel Usner)
    #     # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
    #     table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
    #     cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
    #                                                             timingtable = table_base.format('prob'),
    #                                                             effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_lea_z20.eff.fits'),
    #                                                             timingSigma = 0,
    #                                                             maxRadius=400.,
    #                                                             quantileEpsilon=1e-8,
    #                                                             )
    # elif options.ICEMODEL == 'spice_lea_no_eff_dist':
    #     logging.log_info('Using Spice Lea tables:')
    #     # add table base/cascade service (photonics tables excluding effective distande parametrization by Marcel Usner)
    #     # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
    #     table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_lea_flat_z20_a10.{}.fits')
    #     cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
    #                                                             timingtable = table_base.format('prob'),
    #                                                             timingSigma = 0,
    #                                                             maxRadius=400.,
    #                                                             quantileEpsilon=1e-8,
    #                                                             )
    # elif options.ICEMODEL == 'spice_3.2.1':
    #     logging.log_info('Using Spice 3.2.1 tables:')
    #     # try with the new photonics tables that tianlu produced (end of 2021)
    #     # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
    #     table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_3.2.1_flat_z20_a5.{}.fits')
    #     cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
    #                                                             timingtable = table_base.format('prob'),
    #                                                             effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_3.2.1_z20.eff.fits'),
    #                                                             tiltTableDir = os.path.expandvars('$I3_SRC/photonics-service/resources/tilt/'),
    #                                                             timingSigma = 0,
    #                                                             maxRadius=400.,
    #                                                             quantileEpsilon=1e-8,
    #                                                             )
    # elif options.ICEMODEL == 'spice_3.2.1_no_eff_dist':
    #     logging.log_info('Using Spice 3.2.1 tables:')
    #     # try with the new photonics tables that tianlu produced (end of 2021)
    #     # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
    #     table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_3.2.1_flat_z20_a5.{}.fits')
    #     cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
    #                                                             timingtable = table_base.format('prob'),
    #                                                             tiltTableDir = os.path.expandvars('$I3_SRC/photonics-service/resources/tilt/'),
    #                                                             timingSigma = 0,
    #                                                             maxRadius=400.,
    #                                                             quantileEpsilon=1e-8,
    #                                                             )
    # elif options.ICEMODEL == 'spice_3.2.1_no_tilt':
    #     logging.log_info('Using Spice 3.2.1 tables:')
    #     # try with the new photonics tables that tianlu produced (end of 2021)
    #     # needs this photonics-service project http://code.icecube.wisc.edu/svn/projects/photonics-service/branches/ice-properties
    #     table_base = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_single_spice_3.2.1_flat_z20_a5.{}.fits')
    #     cascade_service = photonics_service.I3PhotoSplineService(amplitudetable = table_base.format('abs'),
    #                                                             timingtable = table_base.format('prob'),
    #                                                             effectivedistancetable = os.path.expandvars('$I3_DATA/photon-tables/splines/cascade_effectivedistance_spice_3.2.1_z20.eff.fits'),
    #                                                             timingSigma = 0,
    #                                                             maxRadius=400.,
    #                                                             quantileEpsilon=1e-8,
    #                                                             )
    else:
        logging.log_info('Selected ice model ({}) currently not supported'.format(options.ICEMODEL))

    excludedDOMs = tray.Add(HighEnergyExclusions,
                            Pulses = pulses,
                            ExcludeDeepCore = None,
                            ExcludeSaturatedDOMs = None,
                            ExcludeBrightDOMs = 'BrightDOMs',
                            BrightDOMThreshold = options.BrightDOMThreshold,
                            SaturationWindows = 'SaturationWindows',
                            BadDomsList = 'BadDomsList',
                            CalibrationErrata = 'CalibrationErrata',
                            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
                            )

    # set TaupedeFit and MuMillipedeFit parameters
    common_params = dict(
        Pulses=pulses,
        CascadePhotonicsService=cascade_service,
        ExcludedDOMs=excludedDOMs,
        minimizer='SIMPLEX',
        )

    logging.log_info('common_params: {}'.format(common_params))

    name = 'TaupedeFit'

    # run the fit (through the TaupedeFitWrapper)
    tray.Add(TaupedeFitWrapper,
            name,
            Seed = Seed,
            # some explicitely set parameters
            Iterations = 1,
            PhotonsPerBin = 1,
            lengthbounds = [-800., 800.],
            If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
            **common_params
            )

    MuMillipede_name = 'MuMillipedeFit'
    MuMillipede_Seed = 'TaupedeFit'

    # run function to correct for negative decay lengths and write reco data as keys
    tray.AddModule(
        correct_negative_decay_lengths,
        'CorrectNegativeDecayLength',
        fitname=MuMillipede_Seed,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    common_params['PhotonsPerBin'] = 1 
    common_params.pop('minimizer')

    # run the fit (with the MuMillipede_Seed (probably TaupedeFit) as Seed)
    tray.AddModule(
        "MuMillipede",
        MuMillipede_name,
        Output=MuMillipede_name,
        ShowerSpacing=5,
        # Smaller shower spacing
        # ShowerSpacing=2,
        MuonSpacing=0,
        SeedTrack = MuMillipede_Seed,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        **common_params
        )

    # add function to write true data as keys
    tray.AddModule(
        write_reco_data_to_keys,
        "WriteRecoDataToKeys",
        fitname=name,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    # add function to write taupede fit parameters
    tray.AddModule(
        write_taupede_fitparameters,
        "WriteTaupedeParamsToKeys",
        fitname=name,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    # add function to extract MuMillipede energies
    tray.AddModule(
        extract_MuMillipede_energies,
        "extractMuMillipedeenergies",
        MuMillipede_fitname=MuMillipede_name,
        Taupede_fitname=name,
        If=lambda frame: frame["I3EventHeader"].sub_event_stream == "InIceSplit",
        )

    if write_hdf5:

        logging.log_info('hdf5 keys to store: {}'.format(HDF5_KEYS))

        outfile_hdf5 = options.OUTFILE.replace('.i3.zst',".hdf5")

        # get weights
        tray.AddModule(
            store_weights,
            "store_weights",
            Streams=[icetray.I3Frame.Physics],
        )

        tray.Add(
            I3HDFWriter,
            Output=outfile_hdf5,
            Keys=HDF5_KEYS,
            SubEventStreams=["InIceSplit"],
        )

    tray.Add("I3Writer",
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        filename=options.OUTFILE,
        )

    tray.Execute()
    tray.Finish()

    t1 = time.time()
    logging.log_info('Time it took: {:.3f}s'.format(t1-t0))
    print('done ...')


if __name__ == '__main__':
    main()
