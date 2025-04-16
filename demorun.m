%%% demo
close all
addpath(genpath(pwd)) 
mex analysis/spktime2count.c 

data_folder = [pwd '/'];

T=3000;  % total time of simulation (ms)

% input = -0.6;
input = 0;
% input = 0.6;

targetpop = 2; % E - 1; PV - 2; SOM - 3; VIP - 4.
filename=strrep(sprintf('%sTestfile_%.2f',data_folder,input),'.',''),

%%%%%%%%%
% tic
SetUp_Sim(filename, data_folder, T, input, targetpop);
% toc

% tic
Sim_Calc(data_folder, filename);
% toc

%tic
Sim_Figure(filename,T, input, targetpop);
%toc

%% functions to collect results across trials and input values 
% resampleCount = 1;
% DataExtraction_CohPwrSpc(data_folder, fnamesave, like_filename, resampcount)
% DataExtraction_FRCorr(data_folder, fnamesave, like_filename)
% DataExtraction_Current(data_folder, fnamesave, like_filename)
% Calc_CohPeakandThreshold(name,data_folder)

