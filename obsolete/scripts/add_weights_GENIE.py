#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.0/icetray-start
#METAPROJECT /data/user/dvannerom/Metaprojects/OscNext/build/

import collections, os, sys, copy, socket
from optparse import OptionParser
from icecube.oscNext.frame_objects.weighting import weighting

from icecube import icetray
from I3Tray import *

from icecube import tableio, hdfwriter

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)
#parser.add_option("-i", "--inputfile", type=str, default='/data/sim/DeepCore/2018/pass2/genie/step1/140000/NuMu_D_140000_000000_step1.i3.gz', help="Input file path")
parser.add_option("-i", "--inputfile", type=str, default='/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001521_level2_converted.zst', help="Input file path")
#parser.add_option("-i", "--inputfile", type=str, default='/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000327_level2_converted.zst', help="Input file path")
parser.add_option("-o", "--outputfile", type=str, default='test.i3.gz', help="Output file path")
(options,args) = parser.parse_args()

#inputfile      = options.inputfile
#inputfile      = ['/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001500_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001501_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001502_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001503_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001504_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001505_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001506_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001507_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001508_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/140000/NuMu_140000_001509_level2_converted.zst'
#                 ]
#inputfile      = ['/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000590_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000591_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000592_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000593_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000594_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000595_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000596_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000597_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000598_level2_converted.zst',
#                  '/data/sim/DeepCore/2018/pass2/genie/level2/120000/NuE_120000_000599_level2_converted.zst'
#                 ]
inputfile      = ['/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000340_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000341_level2_converted.zst', 
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000342_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000343_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000344_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000345_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000346_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000347_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000348_level2_converted.zst',
                  '/data/sim/DeepCore/2018/pass2/genie/level2/160000/NuTau_160000_000349_level2_converted.zst'
                 ]
outputfile     = options.outputfile

tray = I3Tray()

#tray.Add("I3Reader",Filename=inputfile)
tray.Add("I3Reader",FilenameList=inputfile)

# Store the weighting info
tray.Add( weighting, "oscNext_weighting", data_type="genie", dataset=1, overwrite=True, num_files=len(inputfile) )

tray.Add("I3Writer",Filename=outputfile)

tray.Execute()
tray.Finish()

print('done ...')
