function SuppFig15(varargin) 

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = [pwd '/'];

setPlotOpt('custom','path',cDir,'width',13,'height',16); 
outpath=cDir;
Sep = '/';



% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(3,4,[0.1 0.1 0.8 0.85], .35, .35)';% DC = axesDivide(4,1,[0.1 0.2 0.85 0.65], .5, 0.5)';

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
      AH(i) = axes('Pos',DC{i}); hold on; 
end

if P.Recompute  
    LF_generateData(fname);
end
HF_setFigProps;

% START PLOTTING 

figname = 'SuppFig15';

data(1) = load('AllCohPwr_Avgs_resamp_SOM_taud_baseline');
data(2) = load('AllCohPwr_Avgs_resamp_PV_taud_baseline');
data(3) = load('AllCohPwr_Avgs_resamp_PV_taud_subcircuit');

maxyrates = [20 20];

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

Npop = [4 4 2]; 
%%% for powerspectrum
%%%% only plotting a few static input cases in row 3 to avoid overcrowding
%%%% each panel

taucase = {'SOM \tau_d', 'PV \tau_d','PV \tau_d subcircuit'};
XLim =[0 40]; 
for k = 1:3
    %%% Row 1: FR versus static input
    iA=k;
    axes(AH(iA));
    for pop=1:Npop(k)
        plot(data(k).TrialAvgData.param, data(k).CohPeakData.nuSim(:,pop),'color',colororder(pop,:),'linewidth',1.2)
        if k ==1
            text(1.07,1.0-0.15*pop,Pops{pop},'unit','n','color',colororder(pop,:),'FontSize',8);
        end
        if k == 1
            ylabel('Firing rate (Hz)','FontSize',10);
        end
        xlabel('\tau_d (ms)','FontSize',10);
    end
    tt = title(taucase{k},'FontSize',10);
%     xlim(XLim);
%     ylim([0 maxyrates(k)+5]);
    set(gca, 'FontSize',10);
end

for k = 1:3
    %%% Row 2: avg maximum coherence versus static input
    iA=k+3;
    axes(AH(iA));
    for pop=1:Npop(k)
        plot(data(k).TrialAvgData.param, data(k).CohPeakData.firstPeak(pop,:),'color',colororder(pop,:),'linewidth',1.2)
    end
   if k == 1
    ylabel('Max coherence','FontSize',10);
   end
    xlabel('\tau_d','FontSize',10);
%     xlim(XLim);
    ylim([0 1]);
    set(gca, 'FontSize',10);
end

YLim ={[10 20],[10 20],[10 35]};
for k = 1:3
    %%% Row 2: avg maximum coherence versus static input
    iA=k+6;
    axes(AH(iA));
    for pop=Npop(k):-1:1
        plot(data(k).TrialAvgData.param, data(k).CohPeakData.firstPeakFreq(pop,:),'color',colororder(pop,:),'linewidth',1.2)
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


for k = 1:3
  
    %%% Row 4: coherence versus frequency
    iA=k+9;
    axes(AH(iA));
    colors = cool(length(data(k).TrialAvgData.param));

    hold on
    pop = 1; 
    for pid = 1:length(data(k).TrialAvgData.param)
        cohrmean = (data(k).TrialAvgData.Cd_data_trialavg{pid}{pop,pop}...
            *data(k).TrialAvgData.Npair_data_trialavg{pid}{pop,pop})...
            /sum(data(k).TrialAvgData.Npair_data_trialavg{pid}{pop,pop});
        plot(data(k).TrialAvgData.freq, cohrmean, 'Color',colors(pid,:), 'LineWidth',1.15);
        peakfreq = data(k).CohPeakData.firstPeakFreq(pop,pid);
        maxcohr = data(k).CohPeakData.firstPeak(pop,pid);
        plot(peakfreq,maxcohr,'*','color',colors(pid,:))
        text(.8,1.-0.1*pid,num2str(data(k).TrialAvgData.param(pid)),'unit','n','color',colors(pid,:),'FontSize',8);
    end
    xlim([0 50])
    ylim([0 1])
    xlabel('frequency (Hz)')
    ylabel('coherence')
end


HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig']);
set(gcf, 'Renderer', 'painters');
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end

