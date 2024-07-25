#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/

# # metaproject on kepler
# #METAPROJECT /afs/ifh.de/user/l/lfischer/scratch/ICsoft/meta-projects/oscnext_meta/build/


### The Generation Level simulation processing script for simple double cascades ###

#
# Author: Leander Fischer
#

import time
t0 = time.clock()

from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o","--outfile",dest="outfile",type="string", help="Outfile name")
parser.add_option("-s","--seed",dest="seed",type="int", help="Require a unique number for each job/file (use filenumber)")
parser.add_option("-n","--nEvents", dest="nEvents",type="float", help="Number of total events [default: %default]", default = 1000)
parser.add_option("-i","--index",dest="index",type="int", help="Spectral index for energy sampling [default: %default]", default = 2)
parser.add_option("--xcenter", dest="xcenter",type="float", help="Cylinder center x [default: %default]", default = 46.29)
parser.add_option("--ycenter", dest="ycenter",type="float", help="Cylinder center y [default: %default]", default = -34.88)
parser.add_option("--zcenter", dest="zcenter",type="float", help="Cylinder center z [default: %default]", default = -330.0)
parser.add_option("--xyradius", dest="xyradius",type="float", help="Cylinder radius [default: %default]", default = 150.)
parser.add_option("--zheight", dest="zheight",type="float", help="Cylinder height (full) [default: %default]", default = 300)
parser.add_option("--Emin", dest="Emin",type="float", help="Min total energy in GeV [default: %default]", default = 1.)
parser.add_option("--Emax", dest="Emax",type="float", help="Maximum total energy in GeV [default: %default]", default = 10000.)
parser.add_option("--Zmin", dest="Zmin",type="float", help="Min zenith angle in degree [default: %default]", default = 70.)
parser.add_option("--Zmax", dest="Zmax",type="float", help="Maximum zenith angle in degree [default: %default]", default = 180.)
parser.add_option("--Amin", dest="Amin",type="float", help="Min azimuth in degree [default: %default]", default = 0.)
parser.add_option("--Amax", dest="Amax",type="float", help="Maximum azimuth in degree [default: %default]", default = 360.)
parser.add_option("--Lmin", dest="Lmin",type="float", help="Min decay length in meters [default: %default]", default = 0.)
parser.add_option("--Lmax", dest="Lmax",type="float", help="Maximum decay length in meters [default: %default]", default = 1000.)
parser.add_option("--hdf5", action="store_false", dest="write_hdf5", default=True, help="Write hdf5 file [default: %default]")
parser.add_option("-l","--logtofile", action="store_true", dest="logtofile", default=False, help="Write log to file with same name [default: %default]")
options,args = parser.parse_args()


import numpy as np

outfile 	= options.outfile
seed 		= (options.seed)
n_events 	= int(options.nEvents)
gamma       = (options.index)
x_center    = (options.xcenter)
y_center    = (options.ycenter)
z_center    = (options.zcenter)
xy_radius   = (options.xyradius)
z_height    = (options.zheight)
energymin 	= (options.Emin)
energymax 	= (options.Emax)
zenithmin 	= np.radians((options.Zmin))
zenithmax 	= np.radians((options.Zmax))
azimuthmin 	= np.radians((options.Amin))
azimuthmax 	= np.radians((options.Amax))
lengthmin 	= (options.Lmin)
lengthmax 	= (options.Lmax)
write_hdf5  = options.write_hdf5

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

if options.logtofile:
    outfile_log = outfile.replace('.i3.zst',".log")
    assert(outfile_log != outfile)
    handler = logging.FileHandler(outfile_log)
    log.addHandler(handler)
else:
    handler = logging.StreamHandler()
    log.addHandler(handler)


log.info('Outfile: {}'.format(outfile))
log.info('Random seed: {}'.format(seed))
log.info('Events per file: {}'.format(n_events))
log.info('Spectral index: {}'.format(gamma))
log.info('Sampling cylinder center: [{},{},{}]'.format(x_center, y_center, z_center))
log.info('Sampling cylinder radius: {}'.format(xy_radius))
log.info('Sampling cylinder height: {}'.format(z_height))
log.info('Energy range: [{},{}]'.format(energymin, energymax))
log.info('Zenith range: [{},{}]'.format(zenithmin, zenithmax))
log.info('Azimuth range: [{},{}]'.format(azimuthmin, azimuthmax))
log.info('Length range: [{},{}]'.format(lengthmin, lengthmax))
if(write_hdf5):log.info('Also writing hdf5 file')


import math

# import random
# # set random seed once:
# random.seed(seed)


from icecube import dataclasses, icetray, phys_services
from icecube.hdfwriter import I3SimHDFWriter
from icecube.LeptonInjector import direction_vector, direction_angles, sample_power_law, add_cascades_to_tree, random_sign, sample_exp_in_range
from I3Tray import I3Tray

rndserv = phys_services.I3GSLRandomService(seed = seed)

log.info('Starting ...')
########## Functions needed to store as hdf5

HDF5_KEYS = [
    # I3MCTree
    # true first cascade variables
    "casc0_true_x",
    "casc0_true_y",
    "casc0_true_z",
    "casc0_true_energy",
    "casc0_true_zenith",
    "casc0_true_azimuth",
    "casc0_true_time",
    # true second cascade variables
    "casc1_true_x",
    "casc1_true_y",
    "casc1_true_z",
    "casc1_true_energy",
    "casc1_true_zenith",
    "casc1_true_azimuth",
    "casc1_true_time",
    # total true quantities
    'true_decay_length',
    'true_total_energy',
]

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

########## End


##### Define the simple double cascade generator function #####

def simple_double_cascade_generator(
    frame,
    randomservice,
    # sample one position on cylinder
    x_center=x_center,
    y_center=y_center,
    xy_radius=xy_radius,
    z_center=z_center,
    z_height=z_height,
#     sample total energy from power law
    energymin=energymin,
    energymax=energymax,
    gamma=gamma,
#     sample direction uniform
    zenithmin=zenithmin,
    zenithmax=zenithmax,
    azimuthmin=azimuthmin,
    azimuthmax=azimuthmax,
#     sample length exponential
    lengthmin=lengthmin,
    lengthmax=lengthmax,
):
    tree = dataclasses.I3MCTree()
    primary = dataclasses.I3Particle()
    tree.add_primary(primary)


#     sample one position on a cylinder
    radius = randomservice.uniform(-xy_radius, xy_radius)
    angle = randomservice.uniform(0., 2*np.pi)
    
    x_in_cylinder = radius *  np.cos(angle) + x_center
    y_in_cylinder = radius *  np.sin(angle) + y_center
    z_in_cylinder = randomservice.uniform(z_center-z_height/2, z_center+z_height/2)


#     sample decay length
    decay_length = sample_exp_in_range(a=lengthmin, b=lengthmax, u=randomservice.uniform(0., 1.))


#     sample the direction (same for both)
    direction = [
        randomservice.uniform(zenithmin, zenithmax),
        randomservice.uniform(azimuthmin, azimuthmax),
    ]
    directions = [direction, direction]
    this_direction_vector = direction_vector(direction[0], direction[1])


#     sample which cascade is the first/which is the second
    direction_sign = random_sign(randomservice.uniform(0., 1.))

#     calculate connecting vector and position of second cascade
    connecting_vector = direction_sign * decay_length * this_direction_vector
    x_anywhere = x_in_cylinder + connecting_vector[0]
    y_anywhere = y_in_cylinder + connecting_vector[1]
    z_anywhere = z_in_cylinder + connecting_vector[2]


#     add positions to the list in the correct order
    positions = list()
    if direction_sign > 0:
        positions.append([x_in_cylinder, y_in_cylinder, z_in_cylinder])
        positions.append([x_anywhere, y_anywhere, z_anywhere])
    elif direction_sign < 0:
        positions.append([x_anywhere, y_anywhere, z_anywhere])
        positions.append([x_in_cylinder, y_in_cylinder, z_in_cylinder])
    else:
        log.fatal('Something went wrong with the direction sign sampling.')


#     sample total energy from powerlaw
    total_energy = sample_power_law(k_min=energymin, k_max=energymax, gamma=gamma, u=randomservice.uniform(0., 1.))
   

#     split the total energy
    split_energies = randomservice.uniform(0., 1.)
    energies = [total_energy*split_energies, total_energy*(1-split_energies)]


# #     sample both energies from flat distributions (was used for upgoing/horizontal set)
#     energy_0 = sample_power_law(k_min=energymin, k_max=energymax, gamma=gamma, u=randomservice.uniform(0., 1.))
#     energy_1 = sample_power_law(k_min=energymin, k_max=energymax, gamma=gamma, u=randomservice.uniform(0., 1.))
#     energies = [energy_0, energy_1]
#     total_energy = np.sum(energies)


# #     #Note: This was used for up-going/horizontal set
# #     sample the positions on the cylinder (both)
#     positions = list()
    
#     for count in range(2):
#         radius = randomservice.uniform(-xy_radius, xy_radius)
#         angle = randomservice.uniform(0., 2*np.pi)
        
#         x = radius *  np.cos(angle) + x_center
#         y = radius *  np.sin(angle) + y_center
#         z = randomservice.uniform(z_center-z_height/2, z_center+z_height/2)
        
#         positions.append([x, y, z])
                    
#     connecting_vector = np.array(positions[-1]) - np.array(positions[0])
    
# #     calculate decay length
#     decay_length = np.linalg.norm(connecting_vector)

    
#     #Note: This was used for up-going set
# #     flip the positions if they are in the wrong order (depending on direction)
#     if math.acos(np.dot(connecting_vector, this_direction_vector) / decay_length ):
#         positions = [positions[-1], positions[0]]


# #     sample the direction (same for both) (was used for horizontal set)
#     direction = direction_angles(connecting_vector[0], connecting_vector[1], connecting_vector[2])
#     directions = [direction, direction]
        
#     set time of the second particle to distance/c
    decay_time = decay_length/dataclasses.I3Constants.c

    times = [0., decay_time]
    
    frame["I3MCTree"] = add_cascades_to_tree(frame, tree, positions, directions, energies, times)

##### End of generator function definition #####


##### Run the actual tray stuff #####

# instantiate a tray
tray = I3Tray()

tray.Add("I3InfiniteSource", Stream = icetray.I3Frame.DAQ)

tray.Add(
    simple_double_cascade_generator,
    randomservice = rndserv,
    Streams = [icetray.I3Frame.DAQ],
)

# add function to write true data as keys
tray.AddModule(
    write_true_data_to_keys,
    "WriteTrueDataToKeys",
    Streams = [icetray.I3Frame.DAQ],
    )

tray.Add(
    "I3Writer",
    filename = outfile,
    Streams = [icetray.I3Frame.DAQ,]
)

if write_hdf5:

    log.info('hdf5 keys to store: {}'.format(HDF5_KEYS))

    outfile_hdf5 = outfile.replace('.i3.zst',".hdf5")
    assert(outfile_hdf5 != outfile)

    tray.AddSegment(
        I3SimHDFWriter,
        output = outfile_hdf5,
        keys = HDF5_KEYS,
    )

# execute the tray (multiple times)
tray.Execute(n_events)
tray.Finish()
             
# free memory
del tray

##### End of tray stuff #####


log.info('done ...')
t1 = time.clock()
log.info('Time it took: {:.3f} s'.format(t1-t0))
