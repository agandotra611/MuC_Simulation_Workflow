# MuC_Simulation_Workflow

NOTE: iLCSoft environment that can run ddsim and Marlin needed. Can be generated using:
singularity run /cvmfs/unpacked.cern.ch/registry.hub.docker.com/infnpd/mucoll-ilc-framework:1.0-centos8

# Event generation using MADGRAPH 

Generate signal events in MADGRAPH software with: 
1. generate mu+ mu- > vm vm~ h h
2. generate mu+ mu- > mu+ mu- h h
3. generate mu+ mu- > z h h
4. import model sm-full
generate mu+ mu- > h > h h

Generate background events in MADGRAPH:
1. generate mu+ mu- > vm vm~ b b~ b b~
2. generate mu+ mu- > vm vm~ b b~ h
3. generate mu+ mu- > vm vm~ b b~ z

The pythia and run cards are in the respective card directory for each process. The crossx.html file produced has the cross sections for these processes. 


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


Use lctuple_steer.xml to convert reconstructed plots into root files for use. 


# Calculating invariant mass

Use quadCandTree.C to add more relevant variables to tree, including calculated invariant mass. 


# Getting ratio plots 

Use ratio_plot.C to plot the invariant mass in a stack plot of the signal and three backgrounds. 


