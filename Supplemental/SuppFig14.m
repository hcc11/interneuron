function SuppFig14(varargin)

%%%%% Supp Figure 1 compares default circuits with each population receiving static
%%%% input
%%%%% Figure descrition: 2 rows x 4 cols
%%% each col is pop input case (E input, Pv, Som, and Vip)
%%% row 1: (x,y): static input vs FR
%%% row 2: (x,y): static input vs max coh


figname = 'SuppFig14';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',10,'height',8);
%inpath=[cDir ''];
outpath=cDir;
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(3,2,[0.06 0.1 0.85 0.8],0.35,.55)'; %DC = axesDivide(4,1,[0.1 0.2 0.85 0.5], .5, 0.25)';

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
dataCoh(2) = load('AllCohPwr_Avgs_resamp_PVstatic');
dataCoh(1) = load('AllCohPwr_Avgs_resamp_Sigma10_pvstatic');

Iapp = -1:0.1:1;
idx = abs(Iapp)<=1;
XLim = [-1 1];

Sz = 1*(1:21)+5;
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors;
Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};

freq = 0:500;
pplot = [1 4];
for k = 1:2
    iA=pplot(k);
    axes(AH(iA));
    if k == 2
        text(-0.25,1.1,'Spatial Dependent Connections','unit','n','color','k','FontSize',10,'FontWeight','bold')
    end
    if k == 1
        text(-0.25,1.1,'Spatial Independent Connections','unit','n','color','k','FontSize',10,'FontWeight','bold')
    end

    for pop=1:4
        plot(Iapp(idx), dataCoh(k).CohPeakData.nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
        if k ==1
            text(0.8,0.8-0.1*pop,Pops{pop},'unit','n','color',colororder(pop,:),'FontSize',8)
         end
    %    if k == 1
            ylabel('Firing rate (Hz)','FontSize',10)

    %    end
    end
   xlabel('Static Input','FontSize',10)
   xlim(XLim)

    iA = pplot(k) + 1;
    axes(AH(iA)); 
    for pop = 1:4
        plot(Iapp(idx), dataCoh(k).CohPeakData.firstPeak(pop,idx),'color',colororder(pop,:),'linewidth',1.2)
    end
    ylim([0 1]);
    ylabel('Maximum Coherence','FontSize',10)
    xlabel('Static Input','FontSize',10)
   
    iA = pplot(k) + 2;
    axes(AH(iA));

    colors= cool(length(Iapp(idx)));
    for ii = 1:length(Iapp(idx))
        plot(freq, mean(dataCoh(k).TrialAvgData.Cd_data_trialavg{1,ii}{1,1},2),'color',colors(ii,:),'LineWidth',1.05)
    end
    ylabel('Average E Coherence','FontSize',10)
    xlabel('Frequency','FontSize',10)
    xlim([0 50])
    if k == 2
        colormap(colors);
        c = colorbar;
        xx = Iapp(idx);
        c.Ticks = [0 0.5 1.0];
        c.TickLabels = [xx(1) 0 xx(end)];
        c.Label.String = 'Static Input';
        c.Position = [0.93 0.15 0.009 0.2]; %  x poision of bar, y posiiotn of bar, width of bar, length of bar,
    end
end

HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig'])
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);
close all
end

