#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/lfischer/oscnext_hnl_meta/build

#==============================
# D. Vannerom/L. Fischer
#==============================

import time
from os.path import expandvars
import numpy as np
from optparse import OptionParser
import os
from optparse import OptionParser
from os.path import expandvars
import subprocess
import random

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)

parser.add_option("--low",type="int",default = '1',
                  dest="low", help="the low seed value, or smallest file number")
parser.add_option("--high",type="int", default="1000",
                  dest="high", help="the high seed value")
parser.add_option("--outlocation",type="string",
                  dest="outlocation", help="The folder to send the files to")
parser.add_option("--Emin",type="string", default="2",
                  dest="Emin", help="Minimum energy")
parser.add_option("--Emax",type="string", default="10000",
                  dest="Emax", help="Maximum energy")
parser.add_option("--index",type="string", default="2.",
                  dest="index", help="spectral index")
parser.add_option("--nEvents",type="string", default="100000",
                  dest="nEvents", help="Number of events")
(options,args) = parser.parse_args()

Emax 		= str(float(options.Emax))
Emin 		= str(float(options.Emin))
index 		= str(float(options.index))
low 		= options.low
high 		= options.high
nEvents 	= options.nEvents

outlocation = options.outlocation
if not outlocation.endswith('/'):
	outlocation = outlocation + '/'

if not os.path.exists(outlocation):
	os.makedirs(outlocation)
os.chmod(outlocation, 0o775)

print('Min energy: '+ Emin)
print('Max energy: '+ Emax)
print('Spectral index: '+ index)
print('File range: '+ str(low) +'-'+str(high))

mother_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
process_location = os.path.join(mother_path, 'process')
process_gen_file = os.path.join(process_location, 'process_Gen.py')
print(process_gen_file)

for i in range(options.low,options.high+1):
    seed = int(i)
    seed_num = str(seed).zfill(5)
    outfile = os.path.join(outlocation, 'Gen_{}.i3.bz2'.format(seed_num))
    print([process_gen_file,' -s ',str(seed),' -o ',outfile,' --Emax ',Emax,' --Emin ',Emin,' --index ',str(index),' --nEvents ', nEvents])
    subprocess.call([process_gen_file,' -s ',str(seed),' -o ',outfile,' --Emax ',Emax,' --Emin ',Emin,' --index ',str(index),' --nEvents ', nEvents])
