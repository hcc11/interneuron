function SuppFig7(varargin) 

%%%% Supplpemental Fig 7, plots two cases of E input and PV input in a
%%%% E<=>PV subcircuit system analogous to an E-I subcircuit.
%%% Plots the change of FRs for E and PV across static input changes (panel 1)
%%% and coherence measurements for E and PV across static input changes (panel
%%% 2)

figname = 'SuppFig7';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir = [pwd '/'];
% setPlotOpt('custom','path',cDir,'cols',1,'height',17); 
setPlotOpt('custom','path',cDir,'width',8,'height',8); 
outpath=[cDir];
Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(2,2,[0.1 0.1 0.85 0.85], .5, 0.5)';


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
dataCoh(1) = load('AllCohPwr_Avgs_resamp_TwoPop_Estatic');
dataCoh(2) = load('AllCohPwr_Avgs_resamp_TwoPop_Pvstatic');

limitval =1;
Iapp = -1*limitval:0.1:limitval; 
idx = abs(Iapp)<=limitval; 
XLim = [-1*limitval limitval]; 

corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors; 

Pops = {'E', 'PV', 'SOM', 'VIP'};  
Markers = {'o','o','o','o'};  
 
count = 1;
for k = 1:2
    count;
    iA=count;
    axes(AH(iA));
    for pop=1:2
        plot(Iapp(idx), dataCoh(k).CohPeakData.nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2);
        if k ==1
            text(0.1,1-0.15*pop,Pops{pop},'unit','n','color',colororder(pop,:));
        end
    end
    title(['Input to ' Pops{k}]);
    ylabel('Firing rate (Hz)');
    xlabel('static input');
    xlim(XLim);
    count = count+1;

    iA=count;
    axes(AH(iA));
    for pop=1:2
        plot(Iapp(idx), dataCoh(k).CohPeakData.firstPeak(pop,idx),'color',colororder(pop,:),'linewidth',1.2);
    end
    title(['Input to ' Pops{k}]);
    ylabel('Max coherence');
    xlabel('static input');
    xlim(XLim);
    ylim([0 1]);
count = count+1;
end

HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig'])
set(gcf, 'Renderer', 'painters')
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);
close all
end

