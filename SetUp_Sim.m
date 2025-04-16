function [varargout] = SetUp_Sim(filename, data_folder, T, input, targetpop);
cd('simulation/')
%%% mex to compile c code, individual experiment mex calls below.
% mex spktime2count.c


popwithinput = targetpop; % ordered as 1-E, 2-PV, 3-SOM, 4-VIP
alpha = input; %%% external input current

%%% Options for Simulation/Experiements:
%%%% input value is based on input run through the for-loop for ntrials 
%
% Descriptions of the following options are explained near the end 
opt.static = 1; % static input to network, original starting point for project
opt.heteroinput = 0; %%% static input that has variance across the targetted cells
opt.OUinput = 0; %%% OU input generated for popwithinput
opt.savesrr = 0; %%% Poission input to populations
opt.rampinginput = 0; %%% increasing static input to popwithinput in time, run for extra long simulations

    %%%%%%%%%%%%% this part is for running as a job array on cluster %%%%%%%%%%%
    rng('shuffle');
    % AI = getenv('SLURM_ARRAY_TASK_ID'); %% uncommment if running on cluster
    % job_dex = str2num(AI); %% uncommment if running on cluster
    job_dex = 1; %%% comment out if running on cluster
    seed_offset = randi(floor(intmax/10));
    rng(job_dex + seed_offset);

    %%%%%% define parameters %%%%%%%%%%
    % W-matrix to load
    W_fname=sprintf('%sweight_4pop_4',data_folder);

    %%%%%%%%%
    Tburn=500; % burn for inital transient activity
    dt=.05;  % time step (ms)
    Nt=T/dt;
    param(1).ntrials=1; % number of trials of same simulation, marked as ID

    Npop = 4; % number of populations (E, PV, SOM, VIP)
    N_ratio = [4 0.4 0.4 0.2];
    N = 5e4;

    param(1).Nx=2500;
    param(1).Ncell = ceil(N_ratio/sum(N_ratio)* N);
    Nsum = [0 cumsum(param(1).Ncell)];

    %%%% Record neurons and Isyn
    c=0;
    Nrecord = c*ones(1,length(param(1).Ncell));
    param(1).Loc = rand(N,2);  % (x, y) location for each neuron in the recurrent layer
    Nx1=sqrt(param(1).Nx);
    x=(ceil((1:param(1).Nx)/Nx1)-0.5)/Nx1;
    y=(mod((1:param(1).Nx)-1,Nx1)+1 - 0.5)/Nx1;
    param(1).X_Loc = [x(:),y(:)] ;  % (x, y) location for each neuron in the feedforward layer

    %%%% %% %% %% Network Parameters
    % Connection weights
    Jr = [ 30   -90   -120    0; ...
        40  -150    -60    0; ...
        27     0      0  -10; ...
        72     0    -10    0];

    param(1).Jr=Jr;
    param(1).Jx=[120; 250; 0; 0];
    %%%%%%% poisson input, default to zero
    param(1).Jext= [0 0 0 0];% Jext is external input connection strength, used for possion input (i think)
    param(1).rate = [0.0; 0.0; 0.0; 0.0];

    % Synaptic connection probability
    Prr = [ 0.01    0.04    0.03    0; ...
        0.03    0.04    0.03    0; ...
        0.03    0.      0       0.1; ...
        0.01    0.      0.1     0 ];
    param(1).Prr=Prr;
    param(1).Prx=[ .1; .05; 0; 0];
    param(1).Psyn = ones(Npop+1,1);

    % Connection width parameters
    param(1).sigmaRR=[0.1   0.1  0.2    0; ...
        0.1   0.1  0.2    0; ...
        0.2   0    0      0.2; ...
        0.1   0    0.2    0];
    param(1).sigmaRX = .1*ones(Npop,1);

    param(1).taudsyn = [5; 5; 8; 20; 40]; %synatic decay constant % X, E, PV, SOM, VIP
    param(1).taursyn = ones(Npop+1,1); %synatic rise constant

    param(1).Iapp = [0.; 0.; 0.; 0.];
    param(1).Iapp(popwithinput) = alpha;
    paramVal1=alpha; %% this is saved to each file and used in the calculations for what trials to average together across same input

    %%% just for calculations, these are not used
    Jrscale=param(1).Jr/sqrt(N);
    Jxscale=param(1).Jx/sqrt(N);
    wrx=(Jxscale(:,1)).*param(1).Prx*param(1).Nx;
    wrr=(Jrscale).*param(1).Prr.* repmat(param(1).Ncell,Npop,1);

    %%%%%%%%%%%%% simulations %%%%%%%%%%%%%%
        opt.savecurrent=0; % save Isyn
        opt.iaScale = 0; % ramping external input
        opt.poisson_distribution = 1; %% set as default being 1, just generates empty vector

        %%%%% previously set parameters (left unchanged)
        opt.save=1; % save data
        opt.savesx=1;
        opt.saveS1=1; % don't save spike times from Layer 1 (default is 1, save to filename, or s1_fname if specified)
        opt.saveS2=0; % don't save spike times from Layer 2 (default is 1, save to filename)
        opt.CompCorr=0; %# compute correlations
        Nc=500*ones(1,4);  % # of E neurons sampled from Layer 2 & 1
        opt.loadS1=0;
        %     s1_fname=sprintf('%sRF2D3layer_GaborTheta%.04g_sigma_n%.03g_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f',...
        %     data_folder,theta0,sigma_n,Jx(1),Jx(2),inE,0,trial);
        %     s1_fname=strrep(s1_fname,'.','d');
        opt.fixW=1;
        Wseed1=randi(floor(intmax/10));
        Wseed2=randi(floor(intmax/10));
        opt.saveRm=0; % not written for 4 pop
        opt.Layer1only=0;  % default is 0
        opt.saveParam=1;
        opt.plotPopR = 0; 

        p_stim.Nstim=1;
        p_stim.stim_type='Uncorr';
        p_stim.rX=.01; % Rate of neurons in feedforward layer (kHz), size 1xNstim cell of element 1xNsource array
        p_stim.Nsource=1;  % # of sources for global correlation, size 1xNstim

        opt.saveW=0; 
        opt.useWfile=1;

        mexx = 1;

        %%%%%%%%%%%%%%
        %%%%%%%%%%%%%% %%%%%%%%%%%%%%
        %%% parameters to save - leave these
        ParamChange={'filename', filename;'dt', dt; 'T', T; 'Nc',Nc; 'Tburn',Tburn; 'p_stim',p_stim;...
            'N', N; 'Npop', Npop; 'N_ratio', N_ratio; 'wrx', wrx; 'wrr', wrr; 'alpha', alpha};

        ParamChange=[ParamChange;{'param(1).Jr',param(1).Jr; 'param(1).Jx',param(1).Jx;...
            'param(1).Iapp',param(1).Iapp; 'param(1).taudsyn', param(1).taudsyn; 'param(1).sigmaRX',param(1).sigmaRX;...
            'param(1).sigmaRR',param(1).sigmaRR; 'param(1).taursyn', param(1).taursyn;...
            'param(1).Psyn',param(1).Psyn; 'param(1).Prr', param(1).Prr; 'param(1).Prx', param(1).Prx; ...
            'param(1).Nx',param(1).Nx;'param(1).Ncell',param(1).Ncell; 'Nrecord', Nrecord;...
            'param(1).Loc', param(1).Loc; 'param(1).X_Loc', param(1).X_Loc; 'param(1).N', N;...
            'param(1).ntrials', param(1).ntrials; 'param(1).rate', param(1).rate;...
            'param(1).Jext', param(1).Jext; 'param(1).popwithinput', popwithinput}];

        if opt.loadS1
            ParamChange=[ParamChange;{'s1_fname',s1_fname}];
        end

        if opt.useWfile==1
            if exist([W_fname '.mat'], 'file')==0  % if W_fname does not exist, save W in simulation
                opt.saveW=1;
                opt.useWfile =0; 
            end
        end
        if opt.saveW==1||opt.useWfile==1
            ParamChange=[ParamChange;{'W_fname',W_fname}];
            %ParamChange=[ParamChange;{'W_fname',W_fname; 'Wrrfilename', Wrrfilename; 'Wrxfilename', Wrxfilename}];
        end
         if opt.fixW && opt.useWfile==0
            ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
        end
        %%%%%%%%%%%%%% %%%%%%%%%%%%%%

        %%%%%%%%%%%%%% %%%%%%%%%%%%%%
        %%%%%%% based on the choice of experiments/simulations
        %%%%% This is where different values for an experiment must be made
        if opt.static == 1
            if mexx == 1
                mex EIF_MultiPop_Poisson.c
            end
            %%% STATIC INPUT: constant, homogeneous input to one one
            %%% population based on Iapp
        end

        if opt.heteroinput == 1
            if mexx == 1
                mex EIF_MultiPop_heterogIapp.c
            end
            %%% HETEROGEN/QUENCHED INPUT: input that is static but has a
            %%% distributed value across individiual members of the target population
            param(1).rateStDev =[0; 0; 0.2; 0]; %std of input
            input = [0.; 0.; 0.; 0.];
            input(popwithinput) = alpha; % value of input

            param(1).Iapp = zeros(1,N);
            for pop = 1:Npop
                param(1).Iapp(1,Nsum(pop)+1:Nsum(pop+1)) ...
                    = normrnd(input(pop),param(1).rateStDev(pop),...
                    [1,length(Nsum(pop)+1:Nsum(pop+1))]);
            end
            ParamChange =[ParamChange; {'param(1).rateStDev',param(1).rateStDev}];
            ParamChange{15,2} = param(1).Iapp;
        end

        if opt.OUinput == 1
            if mexx == 1
                mex EIF_MultiPop_Iapp_OUinput.c
            end
            %%% OU INPUT: OU input is targetted to popwithinput in c code
            %%% Iapp still acts like external input, set Iapp to 0 if OU is the
            %%% only 'external' input.

            %%%% OU parameters
            param(1).meanOU = 0;%  %mean value
            param(1).sigmaOU = sqrt(0.16*(2*param(1).tauOU));%%%% for variance 0.16 
            param(1).tauOU = 5; % (ms) %time constant for OU, like E

            param(1).Iapp = zeros(1,N);
            param(1).Iapp(1, (Nsum(popwithinput)+1):Nsum(popwithinput+1))=alpha;
            %%% input for Iapp is 1xN so there is flexibility in how to employ it
            %%% across all cells or being defined independent of popwithinput.
            ParamChange =[ParamChange; {'param(1).meanOU', param(1).meanOU;...
                'param(1).sigmaOU', param(1).sigmaOU;...
                'param(1).tauOU', param(1).tauOU}];
            ParamChange{15,2} = param(1).Iapp;
        end

        if opt.savesrr == 1
            if mexx == 1
                mex EIF_MultiPop_Poisson.c
            end
            %%% POISSON INPUT: generates a poisson spk train based on the
            %%% rate and syn strength to each population
            %%% spk train srr measurements is saved as default
            %%% rate and Jext default values are 0's but can be changed
            %%% here below if using poisson input
            ParamChange{30,2} = [0 0 1 0]; %%% FR of poisson spks, param(1).rate
            ParamChange{31,2} = [0 0 0.5 0]; %%% Syn Strength of poisson spks, param(1).Jext
        end

        if opt.rampinginput == 1;
            if mexx == 1
                mex EIF_MultiPop_Poisson_IappofT_Ramping.c
            end
            %%% RAMPING INPUT: this should be run for a a significantly
            %%% longer period of T (~ 100 or 200 sec) as the external input
            %%% changes at timeshift by the rampvalue
            opt.iaScale = 1; % saves iaScale to file (input values in time)
            timeshift = 5000; % time until next ramped input values
            rampvalue = 0.025; % addition to input (make neg to decreaase)

            ParamChange=[ParamChange;{'param(1).timeshift',...
                timeshift; 'param(1).rampvalue', rampvalue}];
        end

        MultiPopSpatialNet_Simulation(opt, ParamChange);
        save(filename,'ParamChange','paramVal1','-append')
        cd('../')
end
