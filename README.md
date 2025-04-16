# interneuron
Codes for the manuscript:   Edwards, M. M., Rubin, J. E., &amp; Huang, C. (2024). State modulation in spatial networks with three interneuron subtypes. bioRxiv.


1. System Requirements:

Combination of C through MATLAB and MATLAB R2021b (9.11) and 2023b, MathWorks.

Simulations run and processed on combination of CNBC Cluster in the University of Pittsburgh and 2021 MacBook Pro 14” with Apple M1 Max (Sonoma 14.1.1)

2. Installation guide

C compiler is gcc, version “Apple clang version 15.0.0 (clang-1500.0.40.1)”

For MATLAB install information and instructions see:
Download and Install MATLAB

Need to add subfolders to path using: 
addpath(genpath(pwd))

And compiled a C file using:  
mex analysis/spktime2count.c 

3. Demo main file: demorun.m

Demo Run instructions - 
demorun.m calls SetUp_Sim.m which contains the simulation parameters, this subsequently calls MultiPopSpatialNet_Simulation.m to simulate the network with targetpop receiving static input.

OUTPUT: Output is a single simulation with the targetpop recieiving static input equal to the value of input selected. The demo calculates the rates, spiking activity, and spike count correlations and produces a figure with 4 subplots (similar to subplots presented in a single column of Figure 2).

INPUT: (note demorun.m run as is, the following input instructions are suggestions for changes for the user). Uncomment single value of input. Select 1,2,3 or 4 for E, PV, SOM, and VIP respectively to receive static input. 

demorun.m simulates the network for 3000 ms. 
After simulation ends, analysis of the network activity begins and is displayed on the output figure.

Expected demo run time (measured on 2021 Apple M1 Max, macOS Sonoma 14.1.1)
            Generating connection weight matrix: ~ 125 seconds (only needed for first simulation)  
	Simulation Time: ~49 seconds
	Calculation Time: ~39 seconds
	Figure Time: ~14 seconds


4. Reproduction Instructions

To reproduce quantitive results and simulations, changes to input type and network parameters must be made in SetUp_Sim.m (i.e. opt.static, opt.OUinput, or param.Jr for changed to network synaptic strength connections, etc).

For initial simulations, a new weight connection matrix (weight_4pop_4.m) must be generated initially and reloaded on subsequent simulations (i.e. initially requires opt.saveW = 1 and opt.useWfile = 0 for first simulation, then for second and following simulations opt.saveW = 0 and opt.useWfile = 1).

Published results are run for T = 15000 ms and averaged across 5 trials. 

To reproduce published figures - run FigureX.m scripts in figure_scripts/:

figure_scripts/ FigureX.m - where X is a number indicating the figure from the paper. .m scripts are ready to run.

Supplemental/  - contains all Supplemental Figure scripts indicated by SuppFigureX.m

simulation/ - contains functions for simulating the network model. 

analysis/ - contains files for analysis of simulation results. 

BE_figTools/ - figure formatting tools, courtesy of Bernhard Englitz. 

animation/ - contains three example state animation files  and .m scripts to generate them.

data/ - data to reproduce published figures. The full data files for generating figures can be downloaded from https://zenodo.org/records/15232542 
