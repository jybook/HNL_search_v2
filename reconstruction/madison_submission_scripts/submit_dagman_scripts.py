#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Author: Leander Fischer
#

import  os
import logging
import subprocess

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


dagman_location = '/scratch/lfischer/taupede/dagman'

for dag_script in os.listdir(dagman_location):
    if dag_script.endswith('.dag'):
        log.info('Submitting dag script: {}'.format(dag_script))

        command_line = 'condor_submit_dag {}'.format(os.path.join(dagman_location, dag_script))
        log.info('Command Line: "{}"'.format(command_line))
        os.system(command_line)  

        # break
