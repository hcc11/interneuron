function Figure8(varargin)

%%%% Figure 8 - External Input to SOM 
%%% Figure 8 description compares different external input types, OU input
%%% and quenched/heterogenous input where the mean of both inputs equals
%%% the static input value (ranging from -1 to 1).
%%%
%%% Different color curves correspond to changes in the varience of the
%%% external inputs.
%%% Note that the grey plot is for the default static input network case
%%% for comparison.
%%%
%%% Figure output - 1 row x 2 cols
%%% panels are E max coh versus E FR
%%% COL 1: OU external input to SOM
%%% COL 2: Quenched external input to SOM


%%%%%   %%%%%    %%%%%   %%%%%
figname = 'Figure8';
cDir = [pwd '/'];

%%%%%%%%%%
%%% load datasets
%%% static input
dataStat(1) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');

%%% quenched/heteogen static input
dataH(1) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_%sstd_heterogeninput','27','02'),'CohPeakData');
dataH(2) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_%sstd_heterogeninput','27','03'),'CohPeakData'); 

dataH(3) = load('AllCohPwr_Avgs_resamp_QuenchedtoSOMvar12','CohPeakData');
dataH(4) = load('AllCohPwr_Avgs_resamp_QuenchedtoSOMvar16','CohPeakData');

%%% OU is static input to SOM, std = 0.2688: external
 dataOUex2(1) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUex_std%s','27','02'),'CohPeakData'); %
 dataOUex2(2) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUex_std%s','27','03'),'CohPeakData'); %

dataOUex2(3) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar12','CohPeakData');
dataOUex2(4) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar16','CohPeakData');

%%%% %% default for comparison
%dataStat
%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',0);  checkField(P,'Recompute',0);
setPlotOpt('custom','path',cDir,'width',9,'height',5);
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
rows = 1; columns = 2;
DC = axesDivide(columns,rows,[0.1 0.2 0.8 0.65], .3, 0.5)';
Labels = {'A','B','C','D','E','F','G','H'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
    AH(i) = axes('Pos',DC{i}); hold on;
    set(0,'DefaultAxesTitleFontWeight','normal');
end
if P.Recompute
    LF_generateData(fname);
end
HF_setFigProps;

%%%%%%%%
Iapp = -1:0.1:1; inputlim = 1.0;
idx = abs(Iapp)<=inputlim;
Sz = linspace(0.5,40,length(Iapp));
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o','o'};
var_label = {'0.04', '0.09', '0.12','0.16'};
dataset_labels = {'Static Input', 'Quenched Input', 'OU input'};
tempcolor = corecolors(3,:)+0.15*[1,1,1];
colrshift = .15;
pop = 1;

axes(AH(2));
colr = [0.6 0.6 0.6]; %%%% default case
plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.7)
data = dataH;
for ss = 1:size(data,2)
    colr = tempcolor-(ss-1)*[1 1 1]*colrshift;
    colr(colr>1) = 1;
    colr(colr<0) = 0;
    scatter(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
    plot(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.5')

end
title(dataset_labels{2});
ylabel('E max coherence');
xlabel(sprintf('%s Firing Rate (Hz)',Pops{pop}));
ylim([0, 0.8]);

clear data;
%%%OU(external), Jse = 27 input
axes(AH(1));
colr = [0.6 0.6 0.6]; %%% default case
plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.7');

data = dataOUex2;
for ss = 1:size(data,2)
    colr = tempcolor-(ss-1)*[1 1 1]*colrshift;
    colr(colr>1) = 1;
    colr(colr<0) = 0;
    scatter(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
    plot(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.5')
    if ss == 1
        text(0.7,0.65,'var','unit','n','color','k');
        text(0.7,0.55,'0','unit','n','color',[0.6 0.6 0.6]);
    end
    text(0.7,0.55-0.1*ss,var_label{ss},'unit','n','color',colr);
end
title(dataset_labels{3});
ylabel('E max coherence');
xlabel(sprintf('%s Firing Rate (Hz)',Pops{pop}));
ylim([0, 0.8]);
clear data

% SAVE FIGURES
set(gcf, 'Renderer', 'painters')
%savefig([cDir figname '.fig']);
HF_viewsave('path',cDir,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);


end


