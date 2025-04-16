function Calc_WframeCurrent(data_folder,filename, wname, samplesize)


%%%%%% Load W frames (ffwd and rec)
%%% W is same for trials if Wframe == 1
%%%% col 1 = postsyn neuron ID
disp('Current Calc : Loading W matrix');
Ncell = [40000 4000 4000 2000];
Nsum = [0 cumsum(Ncell)];
Nx = 2500;

%samplesize = 100;
   
wname = dir([data_folder wname]);
Wframe = load([data_folder wname.name]);
Kout = Wframe.param.Kr;
Kxout = Wframe.param.Kx;
if size(Wframe.Wrr1,2) == 1
   [Wrr1 , Wrf1] = Wremap(Wframe.Wrr1, Wframe.Wrf1, Kout, Kxout, Ncell, Nx);
    save([data_folder wname.name],'Wrr1','Wrf1','-append');
    Wframe = load([data_folder wname.name]);
end

sampids = zeros(1,4*samplesize);
for poptype = 1:4
    sampids(1,1+(poptype-1)*samplesize:poptype*samplesize) = sort(Nsum(poptype)+ randperm(Ncell(poptype),samplesize));
end

% sampids =  [1025, 20337 24467 33440 33808 ...
% 41421 41234 42260 43260 43339 ...
% 45001 45021 45040 45060 45080 ....
% 49965 49966 49978 49980 49998];
% %sampids(1,1:5) = [1234, 14566, 23457, 23456, 38790];

%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%
%if option.singletrialdata == 1
original_files = dir(filename);
file = load([data_folder original_files.name]);
%matlabdata_filename = sprintf("%s%s",data_folder,original_files.name);
matlabdata_filename = filename;

% nametosave = extractBefore(filename,'.mat');


%m = matfile(matlabdata_filename,'Writable',true);
% datafilename = strcat(char(allnames(end)), '_Current');

% % store variable names from loaded file
s1 = file.s1;
sx = file.sx;

Npop = file.ParamChange{8,2};
N = file.param.N;
Nx = file.param.Nx;
Nc = samplesize*ones(1,Npop);
Ncell = file.param.Ncell;
Nsum = [0 cumsum(Ncell)];
nuSim = file.nuSim; % {E, P, S, V}
Loc = file.param.Loc; %%% this records all cell locations in sim

%Ids = sort(file.param.Irecord);

T = file.param.T;
Tburn = file.ParamChange{5,2};
Jr = file.param.Jr;
Jx = file.param.Jx;
taudsyn = file.param.taudsyn;
taursyn = file.param.taursyn;

Iapp = file.param.Iapp;

Prr = file.param.Prr;
Prx = file.param.Prx;


% % % % % % % % % % % % %
% % % Connection IDs by postsyn pop % % %
%samplesize = 100;  %%% col 1 = post syn, col 2 = presyn
% if option.calculations == 1
disp('Sampling output pop neurons');
tic;

[presynE sampidsE] = Wframe_PresynId_known(Wframe.Wrr1, Kout, Ncell, 1, samplesize, sampids(1:1*samplesize));
[presynPV sampidsPV] = Wframe_PresynId_known(Wframe.Wrr1, Kout, Ncell, 2, samplesize, sampids(1+1*samplesize:2*samplesize));
[presynSOM sampidsSOM] = Wframe_PresynId_known(Wframe.Wrr1, Kout, Ncell, 3, samplesize, sampids(1+2*samplesize:3*samplesize));
[presynVIP sampidsVIP] = Wframe_PresynId_known(Wframe.Wrr1, Kout, Ncell, 4, samplesize, sampids(1+3*samplesize:4*samplesize));

presynEx = Wframe_PresynIdX(Wframe.Wrf1, Kxout, Nx, sampidsE); % find presyn IDs connected to (known) sampled postsyn from ffwd layer
presynPVx = Wframe_PresynIdX(Wframe.Wrf1, Kxout, Nx, sampidsPV);

%presynEx(:,2) = presynEx(:,2);
%presynPVx(:,2) = presynPVx(:,2);

%presynE{1,end+1} = presynEx; % append X_presyn outputs to pops with X inputs
%presynPV{1,end+1} = presynPVx;

%%%% reorder sampled ids
presyn_E = reordercell(presynE,sampidsE,samplesize);
presyn_PV = reordercell(presynPV,sampidsPV,samplesize);
presyn_SOM = reordercell(presynSOM,sampidsSOM,samplesize);
presyn_VIP = reordercell(presynVIP,sampidsVIP,samplesize);

presynX_E = reordercell(presynEx,sampidsE,samplesize);
presynX_PV = reordercell(presynPVx,sampidsPV,samplesize);

clear presynPV presynSOM presynVIP presynEx presynPVx Wframe

        InputCurrent_Calc.Ids = [sampidsE, sampidsPV, sampidsSOM, sampidsVIP];
        save(matlabdata_filename,'-append','InputCurrent_Calc');
        toc;

disp('Calculating spk times and current');
% % % % % % %
%%%% calculation of currents

% % % create connectivity weight vector from Input/presyn -> Output/postysyn
E_scaledweights = [Jr(1,:), Jx(1)];
PV_scaledweights = [Jr(2,:), Jx(2)];
SOM_scaledweights = Jr(3,:);
VIP_scaledweights = Jr(4,:);

%presynspks :output is cell, each cell corresponds to presyn cell type {E,
%PV, SOM, VIP, X}. The cell is row vector with sorted spike times

%%%% output of prespks for each sample and current for each sample

%%% grab vector subset - run is memebr on all inputs to get matrix of spks,
%%% sum and then put to single id

taudsynmodx = [taudsyn(2:end); taudsyn(1)]; %% reorder for X to be in last position
tau_r = taursyn(1); % same tau_r for all pops

dt = 0.05;
tinitial = Tburn;
tfinal = T;  % Change if you want to sample/calc current over longer lenth of time (time(end) = tfinal))s
time = tinitial:dt:tfinal; %ms

% %%% previous files
 %tic;
% E_PreSynspks = presynspks(presynE, s1, sx, T, Npop);
%inputcurrent_E_noweight = inputcurrent_convolv(E_PreSynspks, tfinal, tinitial, dt, taursyn(1), taudsynmodx);
% toc;
% 
% tic;
% PV_PreSynspks = presynspks(presynPV, s1, sx, T, Npop);
% inputcurrent_PV_noweight = inputcurrent_convolv(PV_PreSynspks, tfinal, tinitial, dt, taursyn(2), taudsynmodx);
% toc;
% 
% 
% tic;
% SOM_PreSynspks = presynspks(presynSOM, s1, 0, T, Npop);
% inputcurrent_SOM_noweight = inputcurrent_convolv(SOM_PreSynspks, tfinal, tinitial, dt, taursyn(3),taudsyn(1:4));
% toc;
% 
% tic; 
% VIP_PreSynspks = presynspks(presynVIP, s1, 0, T, Npop);
% inputcurrent_VIP_noweight = inputcurrent_convolv(VIP_PreSynspks, tfinal, tinitial, dt, taursyn(4),taudsyn(1:4));
% toc;

%%%% Output info %%%% %%%% %%%%
%%% X_inputcurrent: [samplesize x length(time(1:20:end))] (every 1 ms)
%%% X_EItotalinput: [ 2 x length(time(1:20:end))] (every 1 ms) avg across
%%% all samples; row 1 - sum of excit current, row 2 - sum of inhib current
%%%
%%% Xoutput_Each: [samplesize x 2 ] : each sample  neuron (each row); col1 output FR and
%%% col2 is output avg current
%%%

tic;
[E_inputcurrent E_EItotalinput Eoutput_eachSamp Eoutput_allSamp] = currentandspks_calc(presyn_E, presynX_E, sampidsE, 1, s1, sx, E_scaledweights, taudsynmodx, taursyn, time, dt, samplesize, Iapp(1), Nsum, Nx, [], []);
toc;
tic;
[PV_inputcurrent PV_EItotalinput PVoutput_eachSamp PVoutput_allSamp] = currentandspks_calc(presyn_PV, presynX_PV, sampidsPV, 2, s1, sx, PV_scaledweights, taudsynmodx, taursyn, time, dt, samplesize, Iapp(2), Nsum, Nx, [], []);
toc;
tic;
[SOM_inputcurrent SOM_EItotalinput SOMoutput_eachSamp SOMoutput_allSamp]= currentandspks_calc(presyn_SOM, [], sampidsSOM, 3, s1, [] , SOM_scaledweights, taudsynmodx, taursyn, time, dt, samplesize, Iapp(3), Nsum, Nx, [], []);
toc;
tic;
[VIP_inputcurrent VIP_EItotalinput VIPoutput_eachSamp VIPoutput_allSamp]= currentandspks_calc(presyn_VIP, [], sampidsVIP, 4, s1, [], VIP_scaledweights, taudsynmodx, taursyn, time, dt, samplesize, Iapp(4), Nsum, Nx, [], []);
toc;

%output_eachSamp = {Eoutput_eachSamp, PVoutput_eachSamp, SOMoutput_eachSamp, VIPoutput_eachSamp};

[C_d,COV_d,Cbar,COVbar,daxis,Var,SC,Nc] = corr_d_current(sampids, [E_inputcurrent; PV_inputcurrent; SOM_inputcurrent; VIP_inputcurrent]', Ncell,Loc,Nc,T,Npop);

InputCurrent_Calc.output_AllSamp.output_AllSamp = {Eoutput_allSamp, PVoutput_allSamp, SOMoutput_allSamp, VIPoutput_allSamp};
%InputCurrent_Calc.output_eachSamp = {Eoutput_eachSamp, PVoutput_eachSamp, SOMoutput_eachSamp, VIPoutput_eachSamp};
InputCurrent_Calc.output_Var = [var(Eoutput_eachSamp); var(PVoutput_eachSamp); var(SOMoutput_eachSamp); var(VIPoutput_eachSamp)];
InputCurrent_Calc.output_Mean = [mean(Eoutput_eachSamp,1); mean(PVoutput_eachSamp,1); mean(SOMoutput_eachSamp,1); mean(VIPoutput_eachSamp,1)];
InputCurrent_Calc.input_EICurrent = {E_EItotalinput, PV_EItotalinput, SOM_EItotalinput, VIP_EItotalinput};
InputCurrent_Calc.input_eachSampCurrent = [E_inputcurrent; PV_inputcurrent; SOM_inputcurrent; VIP_inputcurrent];
InputCurrent_Calc.input_C_d = C_d;
InputCurrent_Calc.input_COV_d = COV_d;
InputCurrent_Calc.input_Cbar = Cbar;
InputCurrent_Calc.input_COVbar = COVbar;
InputCurrent_Calc.input_daxis = daxis;
InputCurrent_Calc.input_Var = Var;
InputCurrent_Calc.input_SC = SC; %equal to current_mean
InputCurrent_Calc.input_Nc = Nc;

save(matlabdata_filename,'-append','InputCurrent_Calc');


clear all;

