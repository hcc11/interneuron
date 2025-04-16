function Figure2(varargin)

close all
%%%%% Figure 2: displays different activity in three activity states
%%%% data is from PV static input case, input values of -0.6 (SS), 0
%%%% (baseline and WS), and 0.6 (SA)

%%%%% Figure descrition: 5 rows x 3 cols
%%%% COLS: each col is example state dynamics
%%% col 1: subcircuit (SA)
%%% col 2: weak sync (WS)
%%% col2 3: strong sync (SS)
%
%%%% ROWS: example data
%%% row 1: snapshat/rasters (proportionally sampled, 500 total neurons)
%%% row 2: rate trace
%%% row 3: bar graph of avg FR
%%% row 4: spkcount correlation
%%% row 5: mean pairwise coherence 


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',0);  checkField(P,'View',0);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir =[pwd '/'];

setPlotOpt('custom','path',cDir,'width',16,'height',16);

outpath=[cDir];
Sep = '/';

figname = 'Figure2';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(3,[2 0.6 0.8 1.5 1.5] ,[0.07 0.1 0.85 0.85], .35, [0.2 0.6 0.4 0.5])';


Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
    AH(i) = axes('Pos',DC{i}); hold on;
    set(0,'DefaultAxesTitleFontWeight','normal');
end
if P.Recompute
    LF_generateData(fname);
end

HF_setFigProps;

% START PLOTTING
%%% plotting esthetics 
Sz = 1*(1:21)+5;
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors;
Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};

titles = {'Subcircuit Asynchonous', 'Weakly Synchronous', 'Strongly Synchronous'};
corrd_labels = {'E-E', 'PV-PV', 'SOM-SOM', 'VIP-VIP'};

%%% organize loaded data
cycledata = load('AllFRCorr_Avgs_PVstatic','IndividData'); % static input == -0.6 %% strong sync state
spkcount_corr = load('AllFRCorr_Avgs_PVstatic','CorrFRData_AcrossTrials');
Cohrdata = load('AllCohPwr_Avgs_resamp_PVstatic','TrialAvgData');
pid = [17 11 6]; 


% cases = [17 11 6]; %% static input datacell values for input == +0.6, 0,0. -0.6 %individ data
% cases2 = [17 11 5]; %% trial data
%%%%%%%%%

Iapp = -1:0.1:1; %range of static input
idx = abs(Iapp)<=1;
XLim = [-1 1]; 

freq = 0:500; %range of freq
daxis = 1/80*(1:2:40); %data_cohrIndivid.IndividData.CohPwrSpec_Calc{1,1}.daxis;
Tw = 8; % small window only filter over 850 seconds for plot
t1 = 4250; t2 = 5000; %snapshot of time simulation
newtime = (t1:t2)-200+1;
datatmp = [];
val = zeros(1,3);
for k = 1:3
    re = cycledata.IndividData.reUnsmoothed{1,k};
    re_smoothed=imfilter(re,ones(Tw,1)/Tw); %filter firing rate
    datatmp = [datatmp; max(re_smoothed)];

    val(k) = max(max(re_smoothed(newtime(1):newtime(end),:)),[],2);
end
ymax = max(max(datatmp));
val = 5*floor(val/5);

N = 5e4;
Ncell = [40000 4000 4000 2000];
Nsum = [ 0 cumsum(Ncell)];
binsize = 1;
samplesize = 500;
nullcolor = [0.15, 0.15, 0.15];

count = 1;
for states = 1:3
    % Row 1 : spk rasters
    iA = states;
    axes(AH(iA));
    spkdata = cycledata.IndividData.s1{1,states};
    initialT = 500; finalT = 15000;
   % paramchange = cycledata.IndividData.paramval(states);
    time = initialT:binsize:finalT;

    popcountinit = 1;
    for pop = 1:4
        nPOPcells = Ncell(pop);
        Popspks = spkdata(:,spkdata(2,:) <= Nsum(pop+1) & spkdata(2,:) > Nsum(pop));
        Popspks = Popspks(:, Popspks(1,:)>t1 & Popspks(1,:)<=t2);
        Popratio = (nPOPcells/50000)*samplesize;

        skips = 1;
        num_Popneurons_firing = size(unique(Popspks(2,:)),2); % number of neurons that are firing
        Popneurons_firing = unique(Popspks(2,:)); % Unique neurons IDs that are firing
        popcountfin = popcountinit + Popratio;

        if num_Popneurons_firing>= Popratio 
            % rng(24);
            samplePopspks = randsample(Popneurons_firing,Popratio,false);
            hold on
            for ii = popcountinit:skips:popcountfin-1
                Popsamplespktimes = Popspks(1,Popspks(2,:) == samplePopspks(ii+1-popcountinit));
                scatter(Popsamplespktimes,ones(1,length(Popsamplespktimes))*ii,2,corecolors(pop,:),'filled');
            end

        elseif num_Popneurons_firing < Popratio && num_Popneurons_firing ~=0
            samplePopspks = randsample(Popneurons_firing,num_Popneurons_firing,false);
            hold on
            for ii = popcountinit:skips:length(samplePopspks)
                Popsamplespktimes = Popspks(1,Popspks(2,:) == samplePopspks(ii));
                scatter(Popsamplespktimes,ones(1,length(Popsamplespktimes))*ii,2,corecolors(pop,:),'filled');
            end
        end
        popcountinit = popcountfin;
        clear Popratio Popneurons_firing num_Popneurons_firing Popspks

    end
    % xlabel('Time (ms)','FontSize',11)
    xlim([t1 t2])
    % set(gca, 'XTick',[t1 ceil(t2+t1)/2 t2],'YTick',[],'FontSize',10)
    set(gca, 'XTick',[],'YTick',[])

    if states == 1
        ylabel('Neuron','FontSize',11)
        % text(-0.2,0.25,'Population Rasters','unit','n','Rotation',90,'FontSize',10);
        for pops = 1:4
           text(1.1,0.4+0.1*pops,[Pops{pops}],'unit','n','color',colororder(pops,:), 'FontSize',8);
        end
    end
    text(0.0,1.1,titles{1,states},'unit','n','FontSize',10,'FontWeight','bold');

    %%% Row 2 : rate traces
    iA = 3+states;
    axes(AH(iA));
    re = cycledata.IndividData.reUnsmoothed{1,states};
    re_smoothed=imfilter(re,ones(Tw,1)/Tw); %smooth rate

    hold on
    for pop=1:4
        plot((t1:t2), re_smoothed(newtime(1):newtime(end),pop)','color',colororder(pop,:),'linewidth',1.2)
    end
    if states == 1
        text(-0.2,0.05,'Rate (Hz)','unit','n','Rotation',90,'FontSize',10);
        set(gca, 'XTick', [],'YTick',[0 val(states)])

    end
    xlim([t1 t2])
    xlabel('Time (ms)','FontSize',11)
    set(gca, 'XTick', [t1 ceil(t2+t1)/2 t2],'YTick',[0 val(states)],'FontSize',10)
    hold off
    clear re
   
    % re_each = outputFRcalc(1:50000,15000,250, spkdata);
    Tburn = 250; T = 15000;
    re_each =hist(spkdata(2,spkdata(1,spkdata(2,:)>0)>Tburn),1:5e4)/(T-Tburn)*1e3;
    mean_re = [mean(re_each(1:40000)),mean(re_each(40001:44000)), mean(re_each(44001:48000)), mean(re_each(48001:50000))];
    std_re = [std(re_each(1:40000)),std(re_each(40001:44000)), std(re_each(44001:48000)), std(re_each(48001:50000))];
    stderror_ofmean_re = std_re./sqrt(Ncell);

    %%% Row 3: bar plot for average rates
    iA = 6+states; 
    axes(AH(iA));
    b= bar(cycledata.IndividData.nuSim{1,states},'FaceColor','flat');
    b.CData = corecolors;
    box off
    hold on
    for iii = 1:numel(b.XEndPoints)
        xtips = b.XEndPoints(iii);
        ytips = b.YEndPoints(iii);
        errorbar(xtips,ytips,stderror_ofmean_re(iii),'k.','LineWidth',1.2);
    end
    hold off

    mxval = ceil(max(cycledata.IndividData.nuSim{1,states})/5)*5;
    ylim([0 mxval]);
    xticks([1:4]);
    xticklabels({'E','PV','SOM','VIP'});
    set(gca,'YTick',[0 mxval],'FontSize',10);
    ax = gca;
    ax.XAxis.FontSize = 11;
    if states == 1
        % title('Average Firing Rate','FontSize',12)
        text(-0.2,0.0,'Rate (Hz)','unit','n','Rotation',90,'FontSize',10);
    end
    numsim = round(cycledata.IndividData.nuSim{1,states},2);
    for pops = 1:4
            text(0.9,0.9-0.18*pops,[num2str(numsim(pops),2)],'unit','n','color',colororder(pops,:), 'FontSize',8)
    %        text(1.05,0.9-0.18*pops,['sem:' num2str(round(stderror_ofmean_re(pops),4))],'unit','n','color',colororder(pops,:), 'FontSize',10)
    end

    %%% Row 4: corr_d vs dist
    iA = 9+states; 
    axes(AH(iA));
    hold on;
    for mm = 1:4
        for nn = mm:4
            if size(spkcount_corr.CorrFRData_AcrossTrials.C_d_mean_data_trialavg{1,states}{mm,nn}) == [1 20] & mm == nn
                states
                spkcount_corr.CorrFRData_AcrossTrials.paramVal(states)
                size(spkcount_corr.CorrFRData_AcrossTrials.C_d_mean_data_trialavg{1,states}{mm,nn})
                plot(daxis, spkcount_corr.CorrFRData_AcrossTrials.C_d_mean_data_trialavg{1,states}{mm,nn},'color',corecolors(mm,:),'LineWidth',1.2)

                if states == 1
                    text(0.7,0.9-0.12*mm,[corrd_labels{mm}],'unit','n','color',colororder(mm,:), 'FontSize',8)%, 'FontWeight', 'bold')
                end
            end
        end
    end
    xlabel('Pairwise Distance (a.u.)');
    ylim([-0.05 1]);
    xlim([0 0.5]);
    xticks([0:0.1:0.5]);
 %   ylim([-0.05 1])
 if states == 1
    ylabel('Spike Count Correlation');
 end

% row 5: mean pairwise coherence 
    iA = 12+states; 
    axes(AH(iA));
    hold on;
    for pop = 1:4
        if sum(Cohrdata.TrialAvgData.Npair_data_trialavg{pid(states)}{pop,pop})>1e4
        cohrnorm = Cohrdata.TrialAvgData.Cd_data_trialavg{pid(states)}{pop,pop};%...
            % *Cohrdata.TrialAvgData.Npair_data_trialavg{pid(states)}{pop,pop}(:)...
            % /sum(Cohrdata.TrialAvgData.Npair_data_trialavg{pid(states)}{pop,pop}); 
        [m, I]= max(cohrnorm(2:end)); 
        peakfreq = freq(I+1); 
        maxcohr = m; 
        plot(freq, cohrnorm,'color',colororder(pop,:),'linewidth',1)
        plot(peakfreq,maxcohr,'*','color',colororder(pop,:))
        end
        xlim([0 50])
        xticks([0 25 50])
        xlabel('frequency (Hz)')
        if states == 1
            ylabel('Coherence')
            text(0.7,0.9-0.12*pop,[corrd_labels{pop}],'unit','n','color',colororder(pop,:), 'FontSize',8)%, 'FontWeight', 'bold')
            ylim([0 0.02])
            yticks([0 0.01 0.02])
        elseif states == 2
            ylim([0 0.6])
            yticks([0 0.3 0.6])
        elseif states == 3
            ylim([0 1])
            yticks([0 0.5 1])
        end
    end

    hold off
    count = count +1;
end
HF_setFigProps;

% SAVE FIGURES
set(gcf, 'Renderer', 'painters')
P.Save = 1;
%savefig([outpath figname '.fig'])
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end


