#!/bin/sh /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/icetray-start-standard
#METAPROJECT /data/user/lfischer/software/oscnext_hnl/build/


### The OscNext L3-L5 wrapper processing script ###

'''
Script for running the oscNext event selection.
Enforces:
    CVMFS py environment
    Latest oscNext release (deployed to CVMFS)

Tom Stuttard
'''

if __name__ == "__main__" :

    from icecube.oscNext.selection.oscNext_master import run_oscNext_command_line
    run_oscNext_command_line()
