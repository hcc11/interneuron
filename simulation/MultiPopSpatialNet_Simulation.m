function MultiPopSpatialNet_Simulation(varargin)
% run U & A, with same s1;
%  for stim_type='OriMap_gabor_Tseg'
% save spike counts only
% weights in int32
% rand number generator is 'combRecursive'
% rng(Wseed1,'combRecursive')
% connectivity depends on tuning similarity

% RF2D3layer(option, ParamChange)

% param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr,
%       gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
%       maxns, Irecord, Psyn
%   Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
%   Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii]; % out-degrees are fixed
%   taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
%   sigmaRR=[sigmaee, sigmaei; sigmaie, sigmaii];
%   sigmaRX=[sigmaeX; sigmaiX];

%   Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices,
%       sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
%       cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population.
%       For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and
%       Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop.
%   Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
%       The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.

%   conversion of neuron ID (exc) to (x,y) coordinate in [1, Ne1]x[1, Ne1]:
% exc. ID [1, Ne], x=ceil(I/Ne1); y=(mod((I-1),Ne1)+1); ID=(x-1)*Ne1+y
% inh. ID [Ne+1, Ne+Ni], x=ceil((I-Ne)/Ni1); y=(mod((I-Ne-1),Ni1)+1); ID=(x-1)*Ni1+y+Ne;

% sx: spike trains from Layer0
%     sx(1,:) contains spike times.
%     sx(2,:) contains indices of neurons that spike
% s1: spike trains from Layer1
% s2: spike trains from Layer2
%
% data save in filename
% options is a struct w/ fields:
%   'save','CompCorr','plotPopR','fixW','savecurrent','loadS1','Layer1only'. Default values are 0.
% ParamChange is a cell of 2 columns,
%    the 1st column is variable names and
%    the 2nd column is the values.

% if options.save or options.savecurrent is 1, ParamChange needs to have field 'filename'.
% if options.CompCorr is 1, ParamChange needs to have field 'Nc',e.g. Nc=[500 500];
%      # of neurons to sample from Layer2 & Layer3.
% if options.fixW is 1, ParamChange needs to have field 'Wseed1' & 'Wseed2'.

nVarargs = length(varargin);
switch nVarargs
    case 1
        option = varargin{1};
    case 2
        option = varargin{1};
        ParamChange = varargin{2};
end

if ~isfield(option, 'save') option.save=0; end
%if ~isfield(option, 'CompCorr') option.CompCorr=0; end
if ~isfield(option, 'loadS1') option.loadS1=0; end
if ~isfield(option, 'fixW') option.fixW=0; end
if ~isfield(option, 'useWfile') option.useWfile=0; end
%if ~isfield(option, 'plotPopR') option.plotPopR=0; end
if ~isfield(option, 'savecurrent') option.savecurrent=0; end
if ~isfield(option, 'Layer1only') option.Layer1only=0; end
if ~isfield(option, 'saveRm') option.saveRm=0; end
if ~isfield(option, 'saveSx') option.saveSx=0; end
if ~isfield(option, 'saveS1') option.saveS1=1; end
if ~isfield(option, 'saveS2') option.saveS2=1; end
if ~isfield(option, 'saveParam') option.saveParam=0; end
if ~isfield(option, 'saveW') option.saveW=0; end
%if ~isfield(option, 'plotFIcurve') option.plotFIcurve=0; end
%if ~isfield(option, 'powerSpec') option.powerSpec=0; end
%if ~isfield(option, 'raster1D') option.raster1D=0; end
%if ~isfield(option, 'raster2Dani') option.raster2Dani=0; end
if ~isfield(option, 'savesrr') option.savesrr=0; end
%if ~isfield(option, 'synchronymeasure') option.synchronymeasure=0; end


if option.save==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end
if option.CompCorr==1
    if ~ismember('Nc',ParamChange(:,1))
        error('No Nc (1x2): # of neurons to sample to compute correlations')
    end
end
if option.loadS1==1
    if ~ismember('s1_fname',ParamChange(:,1))
        error('No s1_fname')
    end
end
if option.useWfile==1
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.saveW==0
    if ~ismember('W_fname',ParamChange(:,1))
        error('No W_fname')
    end
end
if option.savecurrent==1
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data')
    end
end

%% define parameters
dim ='2D';

% change parameters
if nVarargs==2
    for i=1:size(ParamChange,1)
        eval([ParamChange{i,1} '= ParamChange{i,2};']);
    end
end

if option.loadS1
    p_stim.s1_fname=s1_fname;
end

% if size(param(2).Iapp,2) ~= size(param(2).Jx,2)
%     error('size(param(2).Iapp,2) should equal size(param(2).Jx,2), which is the number of parameter sets ')
% else
%     Np=size(param(2).Iapp,2);
% end
% for pid=1:Np
%     fprintf('\ninE=%.2f, inI=%.2f\n',param(2).Iapp(1,pid),param(2).Iapp(2,pid))
% end

if length(param(1).sigmaRX) ~= Npop
    error('param(1).sigmaRX should be size Npop x 1')
end
if length(param(1).Jx) ~= Npop
    error('param(1).Jx should be size Npop x 1')
end
%if length(param(1).Iapp) ~= Npop
%    error('param(1).Iapp should be size Npop x 1')
%end
if length(param(1).rate) ~= Npop
    error('param(1).rate should be size Npop x 1')
end

%clear ParamChange varargin;

%% initialization
maxrate=[.05 .035];
param(1).maxns=param(1).N*T*maxrate(1);
fprintf('\nmaximum average rates to record: %d Hz\n',round(maxrate(1)*1e3))

for par=1
    param(par).dt=dt;
    param(par).T=T;
    % EIF neuron paramters
    param(par).gl=[1/15 1/10 1/10 1/10];  % E, PV, SOM, VIP
    %original : param(par).gl=[1/15 1/10 1/10 1/10];  % E, PV, SOM, VIP
    % param(par).gl=[1/15 1/10 ];
    param(par).Cm=ones(1, Npop);
    param(par).vlb= -100 * ones(1, Npop);  % lower bound of Vm
    param(par).vth= -10 *ones(1, Npop); %voltage threshold
    param(par).DeltaT=[2 .5  2  2];
    %      param(par).DeltaT=[2 .5 ];
    param(par).vT= -50 *ones(1, Npop); %mV %voltage threshold
    param(par).vre= -65 *ones(1, Npop); %% reset value
    param(par).tref=[1.5 .5 1.5 1.5];
    %     param(par).tref=[1.5 .5 ];
    V0min=param(par).vre(1);
    V0max=param(par).vT(1);
    param(par).vl = -60 *ones(1, Npop); %leak reversal potential
    rng('shuffle');
    param(par).V0=(V0max-V0min).*rand(param(par).N,1)+V0min;
    
    V0_save = param(par).V0;

    param(par).Kr=ceil(param(par).Prr.* repmat(param(par).Ncell(:), 1, Npop)); % number of postsynaptic cells
    param(par).Kx=ceil(param(par).Prx.* param(par).Ncell');
    Kr = param(par).Kr;
    Kx = param(par).Kx;
    % param(par).Irecord:  index of cells to record
    tmp=[0 cumsum(param(par).Ncell)];
    %param(par).Irecord=[25,250,25500, 30000, 40225, 40176, 40829, 42940, ...
    %    45025, 45673, 45821, 45230, 48134, 48923, 49728, 49283];
    param(par).Irecord =[];
    for pop = 1:Npop
        param(par).Irecord = [param(par).Irecord, randsample((tmp(pop)+1):tmp(pop+1),Nrecord(pop))];
    end
    param(par).Jr=param(par).Jr/sqrt(param(par).N);
    param(par).Jx=param(par).Jx/sqrt(param(par).N);

    %%% Effective connection weights
    % not used in simulation just for calculation
    wrx=(param(par).Jx(:,1)).*param(par).Prx*param(par).Nx;
    wrr=(param(par).Jr).*param(par).Prr.* repmat(param(par).Ncell,Npop,1);
    % For balanced state to exist this vector should be decreasing

%     fprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n\n',wrx(1)/wrx(2),abs(wrr(1,2)/wrr(2,2)),abs(wrr(1,1)/wrr(2,1)));
    % and these values should be >1
%     fprintf('\nAlso, this number should be greater than 1: %.2f\n\n',abs(wrr(2,2)/wrr(1,1)));
    %if par==1
    %    param(par).Imean = p_stim.rX*wrx + param(par).Iapp;
    %    fprintf('\nLayer 1 firing rates for large N: %.2f \n',-(wrr)\(p_stim.rX*wrx + param(par).Iapp)*1e3)
    %end

    if option.OUinput == 1
        param(1).meanOU =  param(1).meanOU;
        param(1).sigmaOU =  param(1).sigmaOU;
        param(1).tauOU = param(1).tauOU;
    end

    if option.rampinginput == 1
        timeshift = param(par).timeshift; % time until next rampped input values
        rampvalue = param(par).rampvalue; % addition to input

        niterations = ceil(param(par).T/timeshift); %number of shifts

        sx_all = cell(1,niterations); % store all data in cells
        s1_all = cell(1,niterations);
        iaScale_all = zeros(Npop,niterations);
        ogIappVal = param(par).Iapp(param(par).popwithinput,1);
        Vm_all = cell(1,niterations);
        Isyn_all = cell(1,niterations);
        Isynprime_all = cell(1,niterations);
        fullT = param(par).T;

        Isyn_initialreload = cell(1,niterations);

        refstate_all = cell(1,niterations);
    end
    param(1).popwithinput = param(1).popwithinput;

end
scurr = rng;

%% generate input spike trains
if option.loadS1==0
    sx=genXspk(p_stim,param(1).Nx,T);
end

%% generate Poisson spike trains in recurrent layer
if option.poisson_distribution == 1 % produce srr poisson input from rr layer
    srr = poisson_distribution(param(1).rate,param(1).Ncell,Npop,T,N);
end

%% Simulation

% Simulate Network
if((param(1).N)<=200000)
    disp('simulation starts')

    save(filename,'T')
    % Random initial membrane potentials
    if option.loadS1
        load(s1_fname,'s1');
        s1=s1(:,s1(2,:)<=param(1).Ne+.1&s1(2,:)>0&s1(1,:)<T);
        disp('load s1')
        %if option.save
        %    save(filename,'T') already done above
        %end
    else
        if option.fixW  %%% W is matrix of connections in reccurent layer (Wrr) and X layer (Wrx)
            if option.useWfile==1
                load(W_fname,'Wrr1','Wrf1')
                data = load(W_fname,'param');
                param.Loc = data.param.Loc;
                fprintf('load weight from %s\n',W_fname)
            else
                param(1).Wseed=Wseed1;
                rng(Wseed1,'twister')
                fprintf('seed%d for Wrr1, Wrf1\n',Wseed1)
                disp('generating Wrr1, Wrf1')
                tic
                [Wrr1,Wrf1, Kout, Kxout]=gen_weights(param(1).Ncell,param(1).Nx,...
                    param(1).sigmaRX,param(1).sigmaRR,param(1).Prr, param(1).Prx,param(1).Loc,param(1).X_Loc);
                elapsetime=toc;
                fprintf('elapsetime=%.2f sec\n',elapsetime)
            end
        else
            disp('generating Wrr1, Wrf1')
            tic
            Wseed1 = rng;
            [Wrr1,Wrf1, Kout, Kxout]=gen_weights(param(1).Ncell,param(1).Nx,...
                param(1).sigmaRX,param(1).sigmaRR,param(1).Prr, param(1).Prx,param(1).Loc,param(1).X_Loc);
            elapsetime=toc;
            fprintf('elapsetime=%.2f sec\n',elapsetime)
        end
        if option.saveW
            save(W_fname,'Wrr1','Wrf1','Wseed1','param', 'Kout', 'Kxout')
            %writematrix(Wrr1, Wrrfilename)
            %writematrix(Wrf1, Wrxfilename)
        end

        disp('simulating Layer1')
        tic

 %%%%% This is where different values for an experiment must be made
 %%
 %%     %%% STATIC INPUT
        if option.static == 1
            [s1,Isyn1,Vm1]=EIF_MultiPop_Poisson(sx, Wrf1,Wrr1,param(1),srr);
        end

 %%     %%% POISSON INPUT
        if option.savesrr==1
            %%% C code called here %%%% %%%%% %%%%%
            [s1,Isyn1,Vm1]=EIF_MultiPop_Poisson(sx, Wrf1,Wrr1,param(1),srr);

            if contains(filename, 'ID1','IgnoreCase',0)==1 & sum(param(1).rate)>0
                save(filename,'srr','-append')
            end
            disp('counting spks srr');
            time = 0:1:T;
            tburn = 500;
            Ncell=  [40000 4000 4000 2000];
            Nsum = [0 cumsum([40000 4000 4000 2000])];
            dt = 0.05; tau_r = 1; tau_d = 5;
            step=1;

            re_srr = zeros(length(time),Npop);
            current_all =cell(1,Npop);

            for pp = find(abs(param(1).rate)>0);%1:Npop
                find(abs(param(1).rate)>0);
                re_srr(:,pp) = hist(srr(1,srr(2,:)<Nsum(pp+1)&srr(2,:)>Nsum(pp)),time)/Ncell(pp)*1e3/step;
                sub_srr = srr(:,srr(2,:)<Nsum(pp+1)&srr(2,:)>Nsum(pp));

                tt = -10*tau_d:dt:10*tau_d; %% filter adapts based on tau_decay
                filter =(exp(-tt./tau_d)-exp(-tt./tau_r)).*(tt>=0);
                filterAll{1,1} = filter/sum(filter);
                count = 1;

                poissoninputIDsamp = sort(Nsum(pp)+ randperm(Ncell(pp),500));
                for nnn = 1:500 %Nsum(pp)+1:Nsum(pp+1)
                    spks_srr_tmp = sub_srr(1,sub_srr(2,:)==poissoninputIDsamp(nnn)); %%%
                    spks_srr_tmp = spks_srr_tmp(spks_srr_tmp<time(end) & spks_srr_tmp>time(1));

                    tmp_current_srr = hist(spks_srr_tmp,time); % ms, size: 1 x T
                    re_tmp_srr = sum(tmp_current_srr)/Ncell(pp); % scalar
                    eachinputrate(count,1) = re_tmp_srr; %rate of each neuron

                    tmp_current_srr = imfilter(tmp_current_srr, filterAll{1,1}, 'conv'); %currrent filter % JJvalue defined above (testvalue = \pm 0.5)
                    %   tmp_current_srr = tmp_current_srr(tburn:step:end).*JJvalue;%scaledweights(prepop);
                    tmp_current_srr = tmp_current_srr(tburn:step:end).*0.5;

                    eachinputcurrent_srr{count,:} = tmp_current_srr; % store current trace

                    avginputcurrent_srr(count,1) = mean(tmp_current_srr); % data to be plotted, comparing srr:F.vs.I and s1: F.vs.I
                    varinputcurrent_srr(count,1) = var(tmp_current_srr);  %% (i.e.) plotting: avginputcurrent(:,col) eachinputrate(:, col)

                    count = count+1; % count so that index storage of data starts at index =1 versus index = nnn = Nsum(pp)+...
                    clear tmp_current_srr spks_srr_tmp re_tmp_srr
                end
                % srr.data.IDsamp{1,1} = poissoninputIDsamp;
                srr_data.eachinputcurrent{1,1}=eachinputcurrent_srr;
                srr_data.avginputcurrent{1,1}=avginputcurrent_srr;
                srr_data.varinputcurrent{1,1}=varinputcurrent_srr;
                clear eachinputcurrent_srr avginputcurrent_srr varinputcurrent_srr
            end
            save(filename,'srr_data','-append')
            %end
            % writematrix(srr(2,:),sprintf('%s_srrID.txt',filename))
            % writematrix(srr(1,:),sprintf('%s_srrTime.txt',filename))
        end


 %%    %%% OU INPUT
        if option.OUinput == 1
            param(1).initialOU=sqrt(param(1).sigmaOU^2/(2*param(1).tauOU))* randn(param(1).Ncell(param(1).popwithinput),1) + param(1).meanOU;
            [s1,Isyn1,Vm1,~,~, ~]=EIF_MultiPop_Iapp_OUinput(sx, Wrf1,Wrr1,param(1),srr);
            
            %%% Fig below requires the recorded last three inputs to check
            %%% OU statistics
%             [s1,Isyn1,Vm1,recordedInput,recordedOU, recordedRandNum]=EIF_MultiPop_Iapp_OUinput(sx, Wrf1,Wrr1,param(1),srr);
%             %
%             fig = figure;
%             subplot(3,2,1)
%             plot(recordedOU);%,Marker='o')
%             xlabel('timestep')
%             ylabel('recorded OU')
%             title(sprintf('mean %.3f, var %.3f',mean(recordedOU),var(recordedOU)))
% 
%             subplot(3,2,2)
%             histogram(recordedOU);
% 
%             subplot(3,2,3)
%             plot(recordedRandNum);%,Marker='o')
%             xlabel('timestep')
%             ylabel('recorded dW')
%             title(sprintf('mean %.3f, var %.3f',mean(recordedRandNum),var(recordedRandNum)))
% 
%             subplot(3,2,4)
%             histogram(recordedRandNum);
% 
%             subplot(3,2,5)
%             plot(recordedInput);%,Marker='o')
%             xlabel('timestep')
%             ylabel('recorded membrane potential ')
%             title(sprintf('mean %.3f, var %.3f',mean(recordedInput),var(recordedInput)))
% 
%             subplot(3,2,6)
%             histogram(recordedInput);

            %   savefig(fig, sprintf('~/Desktop/recorded_OUInput_m%.3f_s%.3f.fig',param(1).meanOU,param(1).sigmaOU))
            % saveas(fig, sprintf('~/Desktop/recorded_OUInput_m%.3f_s%.3f.jpeg',param(1).meanOU,param(1).sigmaOU))
            % close all;

        end

 %%     %%% HETEROGEN/QUENCHED INPUT
        if option.heteroinput == 1
            [s1,Isyn1,Vm1]=EIF_MultiPop_heterogIapp(sx, Wrf1,Wrr1,param(1),srr);
        end

 %%     %%% RAMPING INPUT
        if option.rampinginput == 1
            disp('simulating Layer1: ramping input');
            Nsyn=5;
            param(1).Isyn_init = zeros(1,N*Nsyn); %initialize
            param(1).Isynprime_init = zeros(1,N*Nsyn); %initialize
            param(1).refstate_init= zeros(1,N); %initialize

            s1 = [];
            sx = [];
            Isyn1_rec = [];
            Isynprime1_rec = [];
            Vm1 = [];
            for nn = 1:niterations
                disp(sprintf('iteration =%s',num2str(nn)));

                param(1).Iapp(param(1).popwithinput,1) = (nn-1)*rampvalue + ogIappVal;
                param(1).Iapp(:,1)
                param(1).modnum = timeshift;
                param(1).coeffstep = rampvalue;

                sx_tmp=genXspk(p_stim,param(1).Nx,timeshift);
                param(1).T = timeshift;

                [s1_tmp, Isyn_tmpLASTOUT, Isynprime_tmpLASTOUT, Isyn_tmpFIRSTOUT, Isynprime_tmpFIRSTOUT...
                    Vm_tmp, Vout_tmp, iaScale_tmp,  refstate_tmp, Isyn1_rectmp, Isynprime1_rectmp] ...
                    = EIF_MultiPop_Poisson_IappofT_Ramping(sx_tmp,Wrf1,Wrr1,param(1),srr);

                s1_tmp=s1_tmp(:,1:find(s1_tmp(2,:)==0,1)-1);
                s1_tmp(1,:) = s1_tmp(1,:) + (nn-1)*timeshift;

                sx_tmp(1,:) = sx_tmp(1,:) + (nn-1)*timeshift;

                tmp = zeros(Npop,1);
                unique(iaScale_tmp(1,:));
                unique(iaScale_tmp(2,:));
                unique(iaScale_tmp(3,:));
                unique(iaScale_tmp(4,:));
                tmp(param(1).popwithinput,1) = unique(iaScale_tmp(param(1).popwithinput,:));

                %%%
                s1_all{1,nn} = s1_tmp;
                sx_all{1,nn} = sx_tmp;
                Vm_all{1,nn} = Vm_tmp;

                iaScale_all(:,nn) =  tmp;
                Isyn_all{1,nn} = Isyn_tmpLASTOUT;
                Isynprime_all{1,nn} = Isynprime_tmpLASTOUT;

                Isyn_all{2,nn} = Isyn_tmpFIRSTOUT;
                Isynprime_all{2,nn} = Isynprime_tmpFIRSTOUT;

                Isyn1_rec = cat(2,Isyn1_rec,Isyn1_rectmp);
                Isynprime1_rec = cat(2,Isynprime1_rec,Isynprime1_rectmp);

                Isyn1_recFinal = Isyn1_rectmp(:,end);
                Isynprime1_recFinal = Isynprime1_rectmp(:,end);

                Isyn1_recFirst = Isyn1_rectmp(:,1);
                Isynprime1_recFirst = Isynprime1_rectmp(:,1);

                Vm1 =  Vm_tmp;
                param(1).V0 = Vout_tmp;
                param(1).Isyn_init= Isyn_tmpLASTOUT;
                param(1).Isynprime_init = Isynprime_tmpLASTOUT;
                param(1).refstate_init = refstate_tmp;

                s1 = cat(2,s1,s1_tmp);
                sx = cat(2,sx,sx_tmp);

                clear s1_tmp sx_tmp Vm_tmp Isyn1_rectmp Isynprime1_rectmp Vout_tmp iaScale_tmp tmp  Isynprime_tmpLASTOUT Isyn_tmpLASTOUT
            end
            save(filename, 's1_all', 'sx_all', 'Vm_all', 'Isyn_all', 'Isynprime_all','refstate_all','Isyn1_rec','Isynprime1_rec', '-append' );
        end
        if option.iaScale ==1
            iaScale = iaScale_all;
            save(filename, 'iaScale', '-append')
        end
        if option.savesx ==1
            save(filename, 'sx', '-append')
        end
        if option.savecurrent==1 &&  option.rampinginput == 0
            %%% save postsynaptic current from X, E, PV, SOM, VIP
            %%% postsyn population(receiving pop) is stacked in rows
            %%% Npop*Nrecord = # rows in each Isyn matrix; #cols = T/dt;
            %%% E input current is in 1*Nrecord row block; PV input current is in
            %%% 2*Nrecord row block
            Nskip=1;
            Isyn_X = param(1).Psyn(1)*Isyn1((1:sum(Nrecord)),Nskip:Nskip:end);
            Isyn_E = param(1).Psyn(2)*Isyn1((sum(Nrecord)+1:2*sum(Nrecord)),Nskip:Nskip:end);
            Isyn_PV = param(1).Psyn(3)*Isyn1((2*sum(Nrecord)+1:3*sum(Nrecord)),Nskip:Nskip:end);
            Isyn_SOM = param(1).Psyn(4)*Isyn1((3*sum(Nrecord)+1:4*sum(Nrecord)),Nskip:Nskip:end);
            Isyn_VIP = param(1).Psyn(5)*Isyn1((4*sum(Nrecord)+1:5*sum(Nrecord)),Nskip:Nskip:end);
            Vm1=Vm1(:,Nskip:Nskip:end);

            save(filename,'Isyn_X','Isyn_E','Isyn_PV','Isyn_SOM','Isyn_VIP','Vm1','-append');

            [m,n]=size(Isyn1);
            fprintf('size of Isyn1: [%d, %d]\n',m,n)
            clear Isyn1
        end

        clear Wrr1 Wrf1;
        End=find(s1(2,:)==0,1)-1;
        %         s1=s1(:,s1(2,:)~=0);
        elapsetime=toc;
        clear srr
        tmp=[0 cumsum(param(1).Ncell)];
        for pop= 1:Npop
            nuSim(pop)=1000*nnz(s1(1,1:End)>Tburn & s1(2,1:End)>tmp(pop)& s1(2,1:End)<=tmp(pop+1))/(param(1).Ncell(pop)*(T-Tburn));  % Hz
            fprintf('average rates of pop %d:  %.2f Hz \n',pop, nuSim(pop));
        end
        rate_x=1000*nnz(sx(1,:)>Tburn & sx(1,:)<T)/(param(1).Nx*(T-Tburn));
        fprintf('average rates of input:  %.2f Hz \n', rate_x);
        fprintf('elapsetime=%.2f sec\n',elapsetime)
        wrx=(param(par).Jx(:,1)).*param(par).Prx*param(par).Nx;
        wrr=(param(par).Jr).*param(par).Prr.* repmat(param(par).Ncell,Npop,1);

        if option.save==1
            %save(filename,'sx','-append')
            if strcmp(p_stim.stim_type,'OriMap_gabor_Tseg')
                Tw=50;
                Nt=floor(T/Tw);
                Nt_on=p_stim.T_on/Tw;
                Nt_off=p_stim.T_off/Tw;
                Nseg=Nt_on+Nt_off;
                idx=1:Nt;
                idx(mod(idx-1,Nseg)+1<=Nt_off)=0;
                X=spktime2count(sx,1:param(1).Nx,Tw,Nt,1);
                Xstim=permute(sum(reshape(X(:,idx>.5),param(1).Nx,Nt_on,[]),2),[1 3 2]);
                save(filename,'T','Xstim','Tw')
                if option.saveSx
                    save(filename,'sx','-append')
                end
                clear Xstim X;
                E1=int8(spktime2count(s1(:,1:End),1:param(1).Ne,Tw,Nt,1));
                I1=int8(spktime2count(s1(:,1:End),(1+param(1).Ne):param(1).N,Tw,Nt,1));
                save(filename,'E1','I1', '-append')
                clear E1 I1
            else
                %save(filename,'T')
                %if option.saveSx %not an option
                %    save(filename,'sx','-append')
                %end
            end
            if option.saveRm % needs to be re written for 4 populations model
                %                 re1=hist(s1(1,s1(2,:)<=param(1).Ne&s1(2,:)>0),1:T)/param(1).Ne*1e3;
                %                 ri1=hist(s1(1,s1(2,:)>param(1).Ne),1:T)/param(1).Ni*1e3;
                %                 save(filename,'re1','ri1','-append')
                %                 clear re1 ri1;
            end
            if option.saveS1
                if exist('s1_fname','var')
                    save(s1_fname,'s1')
                    if strcmp(p_stim.stim_type,'OriMap_gabor_Tseg')
                        th_id=p_stim.th_id;
                        save(s1_fname,'th_id','-append')
                    end
                else
                    s1=s1(:,s1(2,:)>0);
                    save(filename,'s1','nuSim','-append')
                end
            end
        end
        s1=s1(:,s1(2,:)<=param(1).N+.1&s1(2,:)>0);
    end

    if option.saveParam
        save(filename,'param','p_stim','wrx','wrr','-append')
    end

else
    error('N too large') % Your computer probably can't handle this
end
disp('simulation ends')

%%
if option.plotPopR
    step=1; 
    time=0:step:T;
    Tw=200;
    Tburn=200;
    re = zeros(length(time),Npop);
    Nsum = [0 cumsum(param.Ncell)];
    for p=1:Npop
        re(:,p)=hist(s1(1,s1(2,:)<=Nsum(p+1)&s1(2,:)>Nsum(p)),time)/param.Ncell(p)*(1e3/step);
    end
    
    re_smoothed=imfilter(re(Tburn/step+1:end,:),ones(round(Tw/step),1)/(round(Tw/step)));
    re_smoothed=re_smoothed(Tw/step/2-1:end-Tw/step/2,:);

    figure
    subplot(2,1,1)
    plot(time,re)
    box off
    xlim([0, T])
    ylabel('FR (Hz)')
    xlabel('time (ms)')

    title('Population rate')
    subplot(2,1,2)
    plot(time(Tburn+Tw/step/2-1:end-Tw/step/2),re_smoothed)
    xlim([0, T])
    box off
    ylabel('FR (Hz)')
    xlabel('time (ms)')
    title('smoothed')
end

