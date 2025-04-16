function SuppFig10(varargin)

%%% Supplemental Figure 10

%%%% Figure description: 4 rows x 6 cols
%%% each plot is the (y,x): E max coherence x E FR
%%% plots the cases of Pv static input and Som static input across modified
%%% circuit strengths
%%% rows are increasing values of Jes = {0, 60, 120, 240}
%%% cols are increasing values of Jps = {0, 60, 120, 240, 360 420}

figname = 'SuppFigure10';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',25,'height',17);
outpath=[cDir];
Sep = '/';

numcols = 7;

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(7,4,[0.03 0.03 0.95 0.9], .5, 0.5)';

indexsubfig_skip = [5, 6, 7, 12, 13, 14, 21];
Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
    if ismember(i,indexsubfig_skip)==1

    else
        AH(i) = axes('Pos',DC{i}); hold on;
        set(0,'DefaultAxesTitleFontWeight','normal');
    end
end
if P.Recompute
    LF_generateData(fname);
end
HF_setFigProps;

% START PLOTTING
% default param Jes = -120, Jps = -60
datafnamesPV={'','AllCohPwr_Avgs_resamp_Jes0_Jps60_pvstatic',...
    'AllCohPwr_Avgs_resamp_Jes0_Jps120_pvstatic','AllCohPwr_Avgs_resamp_Jes0_Jps240_pvstatic_new', '', '', ''; ...
    'AllCohPwr_Avgs_resamp_Jes60_Jps0_pvstatic','AllCohPwr_Avgs_resamp_Jes60_Jps60_pvstatic',...
    'AllCohPwr_Avgs_resamp_Jes60_Jps120_pvstatic','AllCohPwr_Avgs_resamp_Jes60_Jps240_pvstatic', '', '', ''; ...
    'AllCohPwr_Avgs_resamp_Jps0_pvstatic','AllCohPwr_Avgs_resamp_PVstatic',...
    'AllCohPwr_Avgs_resamp_Jes120_Jps120_pvstatic','AllCohPwr_Avgs_resamp_Jes120_Jps240_pvstatic', 'AllCohPwr_Avgs_resamp_Jes120Jps360_pvstatic', 'AllCohPwr_Avgs_resamp_Jes120Jps420_pvstatic',''; ...
    'AllCohPwr_Avgs_resamp_Jes240_Jps0_pvstatic','AllCohPwr_Avgs_resamp_Jes240_Jps60_pvstatic',...
    'AllCohPwr_Avgs_resamp_Jes240_Jps120_pvstatic','AllCohPwr_Avgs_resamp_Jes240_Jps240_pvstatic',  ...
    'AllCohPwr_Avgs_resamp_Jes240Jps360_pvstatic', 'AllCohPwr_Avgs_resamp_Jes240Jps420_pvstatic', 'AllCohPwr_Avgs_resamp_Jes240Jps500_pvstatic'};

datafnamesSOM={'','AllCohPwr_Avgs_resamp_Jes0_Jps60_somstatic',...
    'AllCohPwr_Avgs_resamp_Jes0_Jps120_somstatic','AllCohPwr_Avgs_resamp_Jes0_Jps240_somstatic', '', '', ''; ...
    'AllCohPwr_Avgs_resamp_Jes60_Jps0_somstatic','AllCohPwr_Avgs_resamp_Jes60_Jps60_somstatic',...
    'AllCohPwr_Avgs_resamp_Jes60_Jps120_somstatic','AllCohPwr_Avgs_resamp_Jes60_Jps240_somstatic', '', '', ''; ...
    'AllCohPwr_Avgs_resamp_Jps0_somstatic','AllCohPwr_Avgs_resamp_SOM_static',...
    'AllCohPwr_Avgs_resamp_Jes120_Jps120_somstatic','AllCohPwr_Avgs_resamp_Jes120_Jps240_somstatic', 'AllCohPwr_Avgs_resamp_Jes120Jps360_somstatic', 'AllCohPwr_Avgs_resamp_Jes120Jps420_somstatic', ''; ...
    'AllCohPwr_Avgs_resamp_Jes240_Jps0_somstatic','AllCohPwr_Avgs_resamp_Jes240_Jps60_somstatic',...
    'AllCohPwr_Avgs_resamp_Jes240_Jps120_somstatic','AllCohPwr_Avgs_resamp_Jes240_Jps240_somstatic',...
    'AllCohPwr_Avgs_resamp_Jes240Jps360_somstatic', 'AllCohPwr_Avgs_resamp_Jes240Jps420_somstatic', 'AllCohPwr_Avgs_resamp_Jes240Jps500_somstatic'};

[Jps, Jes] = meshgrid([0 ; -60; -120; -240; -360; -420; -500],[0 ; -60; -120; -240]);

Iapp = -1:0.1:1;
idx = abs(Iapp)<=1;%0.8;
XLim = [-1 1];
Sz = linspace(0.25,40,length(Iapp));
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors;
Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};

pop = 1;
for ii = 1:4 %%% rows
    for jj = 1:numcols %%% cols
        if ii==1&&jj==1
            iA=(ii-1)*numcols+jj;
            axes(AH(iA));
            axis off
            ylim([0 1]);
            xlim([0 1]);

            line([.1 .8], [.1 .1],'Color','k','LineWidth',1.3);
            line([0.2 0.2], [0. .7],'Color','k','LineWidth',1.3);

            text(0.3, -0.1, 'E firing rate (Hz)','FontSize',8, 'FontWeight','bold');
            text(0, 0.15, 'E max coherence','FontSize',8, 'FontWeight','bold','Rotation',90);

            text(0.3, 0.62, 'Input to Pv', 'FontSize',10,'Color',colororder(2,:));
            text(0.3, 0.5, 'Input to SOM', 'FontSize',10,'Color',colororder(3,:));

            count = 1;
        else
            iA=(ii-1)*numcols+jj;
            count = count+1;

            if ismember(count,indexsubfig_skip)==1

            else
                axes(AH(iA));
                data(2) = load(datafnamesPV{ii,jj},'CohPeakData');
                data(3) = load(datafnamesSOM{ii,jj},'CohPeakData');
                for k = 2:3
                    tmp=data(k).CohPeakData.nuSim(:,pop)<0.2;
                    data(k).CohPeakData.firstPeak(pop,tmp)=0;

                    scatter(data(k).CohPeakData.nuSim(idx,pop),data(k).CohPeakData.firstPeak(pop,idx),Sz(idx),colororder(k,:),'filled');
                    plot(data(k).CohPeakData.nuSim(idx,pop),data(k).CohPeakData.firstPeak(pop,idx),'color',colororder(k,:),'linewidth',1);
                end
                title({sprintf('$J_{SOM \\rightarrow E}= $%d',Jes(ii,jj)), sprintf('$J_{SOM \\rightarrow PV}= $%d',Jps(ii,jj))},'Interpreter','latex');
                ylim([0 1]);
            end
        end
    end
end

% SAVE FIGURES
set(gcf, 'Renderer', 'painters');
%savefig([outpath figname '.fig']);
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);
close all
end

