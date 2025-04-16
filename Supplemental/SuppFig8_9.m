function SuppFig8_9(varargin)
close all;


%%% Supplemental Figures 3 and 4:
%%% Figures show the E coherence compared to each population's acriss
%%% synaptic weight changes:
%%%% (suppFig 3) SOM -> VIP; Jvs = -10 (default)
%%%% (suppFig 4) VIP -> SOM; Jsv = -10 (default)
%
%%%% Figure Description - 1 row x 4 cols
%%%% each col is one population in the following order: E, PV, SOM, VIP

%%%%%%%%%%%%%%%%%%

%%%% graphing labels
inhib = {'Jvs', 'SOM to VIP'; 'Jsv', 'VIP to SOM'};
weight_labels = {'0', '-10', '-20', '-50'};


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',0);  checkField(P,'Recompute',0);

cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',16,'height',5);
outpath=[cDir ''];
Sep = '/';


%%%%%%%%%%%%%%%%%%%
%%% load datasets
%%% SOMtoVIP
dataSOM(1) = load('AllCohPwr_Avgs_resamp_SOMtoVIP_Jvs0','CohPeakData');
dataSOM(2) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');
dataSOM(3) = load('AllCohPwr_Avgs_resamp_SOMtoVIP_Jvs20','CohPeakData');
dataSOM(4) = load('AllCohPwr_Avgs_resamp_SOMtoVIP_Jvs50','CohPeakData');

%%% VIPtoSOM
dataVIP(1) = load('AllCohPwr_Avgs_resamp_VIPtoSOM_Jsv0','CohPeakData');
dataVIP(2) = dataSOM(2);
dataVIP(3) = load('AllCohPwr_Avgs_resamp_VIPtoSOM_Jsv20','CohPeakData');
dataVIP(4) = load('AllCohPwr_Avgs_resamp_VIPtoSOM_Jsv50','CohPeakData');

datacases = {dataSOM, dataVIP};
casesFig = {'SOMtoVIP','VIPtoSOM'};

% PREPARE FIGURE
for cc = 1:2
    figname = sprintf('SuppFig%s',num2str(cc+7));
    data = datacases{cc};

    %%%%%%%%
    Iapp = -1:0.1:1; inputlim = 1;
    idx = abs(Iapp)<=inputlim;
    idx(end) = 0;

    Sz = linspace(0.5,40,length(Iapp));
    corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
        0.8500, 0.3250, 0.0980; %p-p
        0.4660 0.6740 0.1880; %s-s %%som is green here
        0.4940, 0.1840, 0.5560]; %v-v

    Pops = {'E', 'PV', 'SOM', 'VIP'};
    Markers = {'o','o','o','o','o'};

    %%%%%%% FR vs coherence
    figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
    rows = 1; columns = 4;
    DC = axesDivide(columns,rows,[0.05 0.2 0.9 0.65], .35, 0.5)';

    Labels = {'A','B','C','D','E','F','G','H'}; LdPos = [-0.1,0.07];
    for i = 1:numel(DC)
        AH(i) = axes('Pos',DC{i}); hold on;
        %      FigLabel(Labels{i},LdPos);
        set(0,'DefaultAxesTitleFontWeight','normal');
    end
    if P.Recompute
        LF_generateData(fname);
    end
    HF_setFigProps;
    %
    colrshift = .2;

    uplim = [7 20 20 20];
    %%%%%%%%
    for pop = 1:4
        tempcolor = corecolors(pop,:)+0.2*[1,1,1];
        sgtitle(inhib{cc,2});
        axes(AH(pop))
        for ss = 1:size(data,2)
            colr = tempcolor-(ss-1)*[1 1 1]*colrshift; colr(colr<0) = 0; colr(colr>1) = 1;
            scatter(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
            plot(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.5');
            if pop == 4
                fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
                if ss == 1
                    text(1.0,1.0-0.12,inhib{cc,1},'unit','n','color','k','FontSize', 10);
                end
                text(1.0,1.0-0.12*(ss+1),weight_labels{ss},'unit','n','color',fontcolr,'FontSize', 8);
            end
        end
        title(Pops{pop});
        ylabel('E Coherence');
        xlabel(sprintf('%s Firing Rate',Pops{pop}));
        xlim([0 uplim(pop)]);

    end

  %  figname = [figname '_FRvCoh'];

    HF_setFigProps;

    % SAVE FIGURES
    set(gcf, 'Renderer', 'painters')
   % savefig([outpath figname '.fig'])
    HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);
    clear figname 
    
    close all
end


end


