function Animation(varagin)
%%% script produces a video of population activity across time (seen as 
%%% single colored does representing a single neuron spiking at that 
%%% moment in time).

%%% Current output runs a 'for-loop' to produce an animation corresponding 
%%% to the three state definitions seen in Fig 2. 

%%% 2D (xy-plane) of spking activity, animation or snapshots at each timepoint
%%% 
addpath(genpath([pwd '/BE_figTools/']));
addpath(genpath([pwd '/data/']));
addpath(genpath([pwd '/animation/']));

cDir = [pwd '/animation/'];
%%%% 
file = load('AllFRCorr_Avgs_PVstatic.mat');

Ncell = [40000 4000 4000 2000]; % num in each population
Loc = file.IndividData.Loc; % cell locations

initT = 4000; %%% ~1000 < initT < finT
              %%% (lower initT bound depends on when first spk happens for recording to start) 
finT = 5000; %%% initT < finT < 15000 (ms)
             %%% (final length simulated time)
statenames = {'SubcircuitAsychnronous', 'WeaklySynchronous', 'StronglySychronous'};
for input = 1:3
    s1 = file.IndividData.s1{1,input}; 
    staticinput = file.CorrFRData_AcrossTrials.paramVal(input);
    savename = sprintf('%s_iT%0.f_fT%0.f',statenames{input},initT,finT);
    
    raster2D_ani(s1,initT,finT,Loc,Ncell,savename,cDir,staticinput);
    clear s1 staticinput savename;
end
end