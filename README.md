# I3_HNL_Decay
Repository for the IceCube HNL+Decay analysis with a low energy "double cascade" signature. The related Wiki can be found [here](https://wiki.icecube.wisc.edu/index.php/Double_Cascade_HNL_Decay).


To build our (custom hnl oscnext) metaproject:

1. Get [this release](https://code.icecube.wisc.edu/svn/sandbox/stuttard/oscNext_meta/releases/V01-00-07) of the OscNext metaproject following the [instructions](https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/stuttard/oscNext/trunk/resources/env/README.md) using the bsm cvmfs `py2-v3.1.1` version `/cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1` (e.g. `source /cvmfs/icecube.opensciencegrid.org/users/BeyondStandardModel/bsm-py2-v3.1.1-1/setup.sh`, source it twice if the error `cp: cannot stat '/etc/OpenCL/vendors/*': No such file or directory` shows up). 

But before creating the build directory do the following steps.

2. Move dataclasses, oscNext, and hdfwriter out of the source folder, e.g.:
```
cd src
mkdir ../replaced_projects
mv dataclasses oscNext hdfwriter ../replaced_projects
```

3. Get our custom version of dataclasses, oscNext, hdfwriter, and LeptonInjector by checking out [this repo](https://github.com/LeanderFischer/LeptonInjector-HNL), e.g. 
```
git init
git remote add origin git@github.com:LeanderFischer/LeptonInjector-HNL.git
git pull origin main
git branch -m main (only needed if your initialised git was named master)
```

4. Now continue to build the software as described in the instructions for the OscNext metaproject using the bsm cvmfs version mentioned above, this is necessary to have the LeptonWeighter module, which is needed for the HNL weighting.

5. After succesfully building you might need to install additional python packages on top of the software build (e.g. tqdm, pyarrow). For this just install any package with `pip install --target=$I3_BUILD/lib package` on top of the already started environment.


The simulation sets are named 19XXXX. 19 is the (PDG like) ID code we have chosen for the sterile neutrino (HNL), 06 is the classifier for our test sets. This is chosen as a unique identifier to distinguish the set names from those used in OscNext (baseline and systematic sets).
