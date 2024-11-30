This is a "dummy" producer that just works (CMSSW_14_0_15)

you can compile it with scram b

you can run it using: cmsRun python/trial_cfg.py

a ROOTFile will be created using a sample few things in the trial_cfg

IMPORTANT
This branch is using alpaka to achieve the very same objective but the class still does not use events from cmssw, rather there is an alpaka kernel that acts as data filler with pseudorandom hits
