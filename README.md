# MuC_Simulation_Workflow

NOTE: iLCSoft environment that can run ddsim and Marlin needed. Can be generated using:
singularity run /cvmfs/unpacked.cern.ch/registry.hub.docker.com/infnpd/mucoll-ilc-framework:1.0-centos8

Gen level hepmc file for mu+ mu- > mu+ mu- h h at /afs/cern.ch/work/a/agandotr/public/3TeV_mumu_mumuhh_pythia8_events.hepmc


# Sim plots

Simulate with command: 

ddsim --steeringFile sim_steer_mumuHbb.py 

Sample simulated file with 100 events at /afs/cern.ch/work/a/agandotr/public/3TeV_mumu_mumuhh_events.slcio

Use LCTuple/examples/lctuple_simhits.xml to generate plots using 

Marlin lctuple_simhits.xml


# Reco plots 

Use the file reco_steer.xml with the 3TeV_mumu_mumuhh_events.slcio as input to perform reconstruction without any beam induced background. Use: 

Marlin reco_steer.xml 

This should produce an output file 3TeV_mumu_mumuhh_RECO.slcio 

To verify whether the reconstruction has been completed, one can run 

anajob 3TeV_mumu_mumuhh_RECO.slcio

which gives event by event information of the slcio file - and one should ideally see non-zero reconstructed particles in each event. 
