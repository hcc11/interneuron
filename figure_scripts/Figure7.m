function Figure7(varargin) 
% Figure 7: Impacts of the spatial and temporal scales of SOM connections on network synchrony.

set(0,'DefaultAxesTitleFontWeight','normal');

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
figname = 'Figure7';
cDir = [pwd '/'];

setPlotOpt('custom','path',cDir,'width',17,'height',16); 
% inpath=[pwd '/data/'];
% inpath1=[cDir 'code/Demo/data/'];

outpath=cDir;

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,3,[0.1 0.1 0.85 0.75], .35, .5)';% DC = axesDivide(4,1,[0.1 0.2 0.85 0.65], .5, 0.5)';
DC2 = axesDivide(3,3,[0.1 0.1 0.85 0.75], .5, .5)';% DC = axesDivide(4,1,[0.1 0.2 0.85 0.65], .5, 0.5)';

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:4
      AH(i) = axes('Pos',DC{i}); hold on;  
end
for i = 5:10
      AH(i) = axes('Pos',DC2{i-1}); hold on;  
end

if P.Recompute  
    LF_generateData(fname);
end
HF_setFigProps;

% START PLOTTING 
plottingpop = 2; % input to PV 
cases = {'default', 'sigmaSE', 'OtherSigmaSOM','AllsigmaSOM'}; 
% figname = 'Figure7';
sigma = 0.05; 
% figname = 'Figure_sigmacases0d1_pvinput';
% sigma = 0.1; 
% sigma_case = 'sigmaSE'; % projection width E->S
% sigma_case = 'OtherSigmaSOM'; %  VIP->S,S->E, S->PV, S->VIP
% sigma_case = 'sigmaFromSOM'; % projection widths from SOM (S->E, S->PV, S->VIP)
% sigma_case = 'AllsigmaSOM'; % all projection widths to and from SOM (E->S, VIP->S,S->E, S->PV, S->VIP) 
if sigma == 0.05
    data(1) =  load('AllCohPwr_Avgs_resamp_PVstatic','CohPeakData','TrialAvgData');
    data(2) = load('AllCohPwr_Avgs_resamp_SigmaSE_0d05_PVinput','CohPeakData','TrialAvgData');
    data(3) = load('AllCohPwr_Avgs_resamp_OtherSigmaSOM_0d05_PVinput','CohPeakData','TrialAvgData');
    data(4) = load('AllCohPwr_Avgs_resamp_narrowSigma2_PVstatic','CohPeakData','TrialAvgData');
elseif sigma == 0.1
    data(1) =  load('AllCohPwr_Avgs_resamp_PVstatic','CohPeakData','TrialAvgData');
    data(2) = load('AllCohPwr_Avgs_resamp_SigmaSE_0d1_PVinput','CohPeakData','TrialAvgData');
    data(3) = load('AllCohPwr_Avgs_resamp_OtherSigmaSOM_0d1_PVinput','CohPeakData','TrialAvgData');
    data(4) = load('AllCohPwr_Avgs_resamp_equalSigma2_PVstatic','CohPeakData','TrialAvgData');
end


limitvalue = 1;
Iapp = -1:0.1:1; 
idx = abs(Iapp)<=limitvalue; 
XLim = [-1*limitvalue limitvalue]; 

%%% plotting esthetics 
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors; 
Pops = {'E', 'PV', 'SOM', 'VIP'};  


%%% for powerspectrum
%%%% only plotting a few static input cases in row 3 to avoid overcrowding
%%%% each panel
freq = 0:500; 
% cohplotcases = {[1 4 8 12 15],[1 4 8 12 15],[1 4 8 12 15],[1 4 8 12 15]};
cohplotcases = {1:2:21,1:2:21,1:2:21,1:2:21, 1:2:21};
% caseshere = sort(unique([1:2:21, 8]));
% cohplotcases = {caseshere, caseshere, caseshere, caseshere};

% taucase = {': default sigma', ': sigma=0.1', ': sigma=0.05'};
markers = {'.','o','d','s'};
for k = 1:4
    %%% Row 1: FR versus static input
    iA=k;
    axes(AH(iA));
    for pop=1:4
        plot(Iapp(idx), data(k).CohPeakData.firstPeak(pop,idx),'color',colororder(pop,:),'linewidth',1.2)
        if k ==1
            text(.9,1.0-0.15*pop,Pops{pop},'unit','n','color',colororder(pop,:),'FontSize',8);
        end
    end
   if k == 1
    ylabel('Max coherence','FontSize',10);
   end
    xlabel('Static Input','FontSize',10);
    xlim(XLim);
    ylim([0 1]);
    set(gca, 'FontSize',10);
    
    tt = title(['Input to ' Pops{plottingpop} ': ' cases{k}],'FontSize',10);
end

data_tau(1) = load('AllCohPwr_Avgs_resamp_SOM_taud_baseline');
data_tau(2) = load('AllCohPwr_Avgs_resamp_PV_taud_baseline');
data_tau(3) = load('AllCohPwr_Avgs_resamp_PV_taud_subcircuit');
Npop = [4 4 2]; 
taucase = {'SOM \tau_d', 'PV \tau_d','PV \tau_d subcircuit'};

for k = 1:3
    %%% Row 2: avg maximum coherence versus static input
    iA=k+4;
    axes(AH(iA));
    for pop=1:Npop(k)
        plot(data_tau(k).TrialAvgData.param, data_tau(k).CohPeakData.firstPeak(pop,:),'color',colororder(pop,:),'linewidth',1.2)
    end
   if k == 1
    ylabel('Max coherence','FontSize',10);
   end
    xlabel('\tau_d','FontSize',10);
%     xlim(XLim);
    ylim([0 1]);
    set(gca, 'FontSize',10);
    tt = title(taucase{k},'FontSize',10);
end

YLim ={[10 20],[10 20],[10 35]};
for k = 1:3
    %%% Row 2: avg maximum coherence versus static input
    iA=k+7;
    axes(AH(iA));
    for pop=Npop(k):-1:1
        plot(data_tau(k).TrialAvgData.param, data_tau(k).CohPeakData.firstPeakFreq(pop,:),'color',colororder(pop,:),'linewidth',1.2)
    end
   if k == 1
    ylabel('Max coherence freq.','FontSize',10);
   end
    xlabel('\tau_d','FontSize',10);
%     xlim(XLim);
%     ylim([0 1]);
    set(gca, 'FontSize',10);
    ylim(YLim{k});
end



HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig']);
set(gcf, 'Renderer', 'painters');
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);


end

