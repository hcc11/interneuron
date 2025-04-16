function Sim_Calc(data_folder, filename)


samplesize = 500; %number of neurons sampled in each pop

wname = 'weight_4pop_4.mat'; % change name if needed

resampleCount = 1; %
%%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%

%%% rate calculations
Calc_corrFR(data_folder, [filename '.mat']) 
%%% coherence calculations
Calc_cspd_d(data_folder, [filename '.mat'], samplesize, resampleCount);

%%% current calculation
%%% - current calculation takes a decent amount of time, best to
%%% calcuate seperatly and first run FR and CSPD calcs for initial
%%% statistics
% Calc_WframeCurrent(data_folder,[filename '.mat'], wname, samplesize)

%%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%  %%%%%%%%%%%%