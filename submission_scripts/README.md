## Instructions to run HNL MC generation&processing

It's possible to run the processing on NPX and for some levels (Gen, Phot) it's also possible on the OSG. Follow the specific instructions to run on NPX/OSG.


### General information for submission

For both NPX and OSG a couple of folders have to be created on the submit node. e.g. `@submit-1` (NPX) or `@pyglidein1-submit`(OSG). On those nodes create a folder in your scratch where you want the `.dag` files to be located (e.g. `/scratch/username/dagman`) and some folders too write log/error/output files (e.g. `/scratch/username/[condor_logs, condor_error, condor_output]`).
This **is different** than the `/scratch/username` on the cobalts, so once again, make sure to be on the submit node to create them. This only needs to be done once and, eventhough the files will be cleared after some time, but the folder structure stays.

**NPX:** To produce the `.dag` files and submit them you **need to be on a submit node** (e.g. `ssh submitter` from the cobalts). Make sure that it shows you being on `@submit-1`.

**OSG:** To produce the `.dag` files and copy all needed files over to the submit node (to submit them from there later) you **need to be on the cobalts** (e.g. `\data\user`) and only later move to the submit node with `ssh sub-1`. To create and copy the files over you will need to also create the `/scratch/username/dagman` directory on cobalt.


### Submit scripts

For every processing script there is a submit script piping the arguments into the processing script.
Unfortunately, the location of the processing script and the location of the log/error/output folders are hardcoded here, so they need to be changed to your respective paths.
Take a look at [submit-Gen.condor](https://github.com/LeanderFischer/I3_HNL_Decay/blob/master/submission_scripts/submit/submit-Gen.condor) for example.
This needs to be done for all the `.condor` scripts in the [submit folder](https://github.com/LeanderFischer/I3_HNL_Decay/blob/master/submission_scripts/submit/), there are specific scripts for the OSG, but changing all is probably a good idea.


### Job scripts

The scripts to create the dag files are located in [jobs](https://github.com/LeanderFischer/I3_HNL_Decay/tree/master/submission_scripts/jobs).
Navigate there on the submit node (for NPX) or on cobalt (for OSG).
The order at which the processing has to be run is the following:

**Gen -> Phot -> Det/L1/l2 -> L3/L4/L5 -> FLERCNN(L6/L7/L8)**

You can take a look at the possible options of each file with (for example) `./submit_Gen.py -h`.
All the defaults are set to reasonable values, but some **have** to be set. For the generation level it is the identifier (new set name) and the HNL mass.
For the higher levels it's just the identifier (to specify which folder to read from and write to) and the `--REALRUN` flag, to write the dag for all files. To run on the OSG you will have to set your username as well as the flag for the OSG and dagman folder (location path should match on cobalt and the submit node).

Once you have chosen all the arguments  you can run the script.

**OSG:** For the OSG you will make sure the copying worked and then go to the submit node and do the following there.

You can take a quick look at the dag file `less /location/of/current_dag_file.dag` to make sure it makes sense.
Additionally, the folder/subfolders to write the simulation i3/hdf5 files should have been created in `/data/ana/BSM/HNL/MC/XXXXXX` (also a good thing to check everything is ok). A metadatada text file with the processing settings should also be found there.

**OSG:** Before submitting to the OSG you have to spin off a grid proxy with `grid-proxy-init -bits 1024 -valid 168:00:00`, where the last argument says how long it should be valid (only works if you have set up your grid certificate following [these instructions](https://wiki.icecube.wisc.edu/index.php/OSG_certs)). It's advised to check occasionally that your proxy is still valid with `voms-proxy-info`. 

If all looks fine you can start the jobs with `condor_submit_dag /location/of/current_dag_file.dag`.
With `condor_q` you can check how the processing is going, every batch (for every submitted dag) shows up as one line with number of done/run/idle/total(/held) jobs.



There is some metadata created by condor with every dag script that is submitted.
Especially `/location/of/current_dag_file.dag.dagman.out` can be helpful in case jobs get held or so.

If some jobs fail and the dag script has finished a rescue dag is created and with `condor_submit_dag /location/of/current_dag_file.dag` (exactly the same command) only the failed jobs will be restarted, automatically. This also works after process/submit scripts have been edited, in case there was some issue.
