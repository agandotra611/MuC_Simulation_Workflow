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

Upon generating each signal/background, store the output in a new directory and begin simulation using 
output name_of_directory/
launch name_of_directory/

This will give a menu where pythia8 shower will be turned OFF. Enter 1 to turn on the pythia shower, and 0 to proceed. 

One can then proceed to edit the param_card and run_card by entering 1 or 2 respectively. The pythia and run cards are in the respective card directory for each process. The crossx.html file produced has the cross sections for these processes. 

This will produce .hepmc files with the processes you need. 


# Sim plots

The hepmc files from MADGRAPH are input for sim_steer_mumuHbb.py 

Simulate with command: 

ddsim --steeringFile sim_steer_mumuHbb.py 

Use LCTuple/examples/lctuple_simhits.xml to generate plots using 

Marlin lctuple_simhits.xml


# Reco plots 

Use the file reco_steer.xml with the .slcio output of sim_steer_mumuHbb.py as input to perform reconstruction without any beam induced background. Use: 

Marlin reco_steer.xml 

This should produce an output file 3TeV_mumu_mumuhh_RECO.slcio 

To verify whether the reconstruction has been completed, one can run 

anajob 3TeV_mumu_mumuhh_RECO.slcio

which gives event by event information of the slcio file - and one should ideally see non-zero reconstructed particles in each event.

Use lctuple_steer.xml to convert reconstructed plots into root files for use. 


# Calculating invariant mass

Use quadCandTree.C to add more relevant variables to tree, including calculated invariant mass. Use:
root -l
.L quadCandTree.C
.x quadCandTree.C("file.root")

# Getting ratio plots 

Use ratio_plot.C to plot the invariant mass in a stack plot of the signal and three backgrounds. 

root -l
.L ratio_plot.C
overlay()


