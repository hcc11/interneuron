function Figure4(varargin) 

%%%% Figure 4 displays changes in one populations FR in all four of the input
%%%% cases compared to the E population coherence

%%%%% Figure description: 1 row x 4 columns
%%%% COLS: each population's FR compared with E Max Coherence for all 4
%%%% cases (i.e. a single panel contains one populations changes in FR from
%%%% each of th four input cases
%%% note - colors of curves correspond to input cases
%%% col 1: E FR verses E max coherence
%%% col 2: PV FR verses E max coherence
%%% col 3: SOM FR verses E max coherence
%%% col 4: VIP FR versus E max coherence

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = [pwd '/'];

setPlotOpt('custom','path',cDir,'width',16,'height',5); 
outpath=[cDir];
Sep = '/';

figname = 'Figure4';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,1,[0.06 0.2 0.85 0.6], [.85 0.5 0.5 -0.4], 0.5)';

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
%       FigLabel(Labels{i},LdPos); 
set(0,'DefaultAxesTitleFontWeight','normal');
end
if P.Recompute  
    LF_generateData(fname);
end
HF_setFigProps;

% START PLOTTING 
data(1) = load('AllCohPwr_Avgs_resamp_Estatic','CohPeakData');
data(2) = load('AllCohPwr_Avgs_resamp_PVstatic','CohPeakData');
data(3) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');
data(4) = load('AllCohPwr_Avgs_resamp_VIPstatic','CohPeakData');

limitval =1;
Iapp = -1*limitval:0.1:limitval; 
idx = abs(Iapp)<=1; 
xlimXLim = [-1*limitval limitval]; 

%%% plotting esthetics 
Sz = linspace(0.25,40,length(Iapp));
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors; 
Pops = {'E', 'PV', 'SOM', 'VIP'};  
Markers = {'o','o','o','o'};  
poplimit=[8, 20, 30, 62]; % to set axes in panel
for pop = 1:4
    iA=pop;
    axes(AH(iA));    
    k = 2;
        scatter(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),Sz(idx(2:end)),colororder(k,:),'filled')
        plot(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),'color',colororder(k,:),'linewidth',1.1)
    k = 4;
        scatter(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),Sz(idx(2:end)),colororder(k,:),'filled')
        plot(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),'color',colororder(k,:),'linewidth',1.1)

    k = 3;
        scatter(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),Sz(idx(2:end)),colororder(k,:),'filled')
        plot(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),'color',colororder(k,:),'linewidth',1.1)
    k = 1;
        scatter(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),Sz(idx(2:end)),colororder(k,:),'filled')
        plot(data(k).CohPeakData.nuSim(idx(2:end),pop),data(k).CohPeakData.firstPeak(1,idx(2:end)),'color',colororder(k,:),'linewidth',1.1)
    if pop == 1
        for k = 1:4
            text(0.84,0.75-0.1*k,['Input to ' Pops{k}],'unit','n','color',colororder(k,:),'FontSize',8);
        end
    end
    xlim([0 poplimit(pop)]);
    
    title([Pops{pop} ' population']);
    xlabel([Pops{pop} ' Firing rate (Hz)']);
    ylabel('E Max coherence');
end 

HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig'])
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end

