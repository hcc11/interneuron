function SuppFigs12_16_17(varargin)

%%%% Supplemental Figures 12, 16, and 17

%   %   %
%%% RECURRENT INPUT FIGS (Jse modifications):
% Fig 12: OU to match I_(E->SOM) as input to SOM and static input to SOM
%          - In single panel, different curves correspond to changes in Jse
% 
%   %   %

%%% EXTERNAL INPUT FIGS:
% Fig 16: OU external input to SOM 
%          - In single panel, different curves correspond to changes in Var
%          of OU input
%
% Fig 17: Quenched external input to SOM 
%          - In single panel, different curves correspond to changes in Var
%          of quenched input
% 
%   %   %

%%%% General Figure Descriptions: 2 rows x 4 cols
%%% COLS - each population's data ordered as E, PV, SOM, VIP
%%% ROWS - 
%%% row 1 - FR versus external input to SOM
%%% row 2 - Max E coh versus FR
%%%
%%% Default circuit is plotted in solid grey line without input markers 

figname = 'SuppFigure'; 
external = 0; %%% Supplemental figures 16 and 17
recurrent = 1; %%% Supplemental figure 12

%%%%%%%
%%% load datasets
%%%% static input
dataStat(1) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData'); 
dataStat(2) = load('AllCohPwr_Avgs_resamp_Jse14_somstatic','CohPeakData');
dataStat(3) = load('AllCohPwr_Avgs_resamp_Jse22_somstatic','CohPeakData');

%%%% heteogenous static input
dataH(1) = load('AllCohPwr_Avgs_resamp_SOM_Jse27_02std_heterogeninput','CohPeakData');
dataH(2) = load('AllCohPwr_Avgs_resamp_SOM_Jse27_03std_heterogeninput','CohPeakData');
dataH(3) = load('AllCohPwr_Avgs_resamp_QuenchedtoSOMvar12','CohPeakData');
dataH(4) = load('AllCohPwr_Avgs_resamp_QuenchedtoSOMvar16','CohPeakData');


%%%% OU external input ===> time and target varying input
%%%% OU input: mean == static input, var == std of E->S current 
dataOUexStat_27(1) = load('AllCohPwr_Avgs_resamp_SOM_Jse27_OUex_std02','CohPeakData');
dataOUexStat_27(2) = load('AllCohPwr_Avgs_resamp_SOM_Jse27_OUex_std03','CohPeakData');
dataOUexStat_27(3) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar12','CohPeakData');
dataOUexStat_27(4) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar16','CohPeakData');


%%%% OU recurrent input ===> like E input current and + static input to SOM
%%%% OU input: mean == mean of E->S current (= 0.65) , var == std of E->S current 
dataOUrecE(1) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','0'),'CohPeakData');
dataOUrecE(2) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','14'),'CohPeakData');
dataOUrecE(3) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar12_Jse20','CohPeakData');
dataOUrecE(4) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','27'),'CohPeakData');

%%%%%%%%%%%%%%%%%%
%%%% graphing labels
var_label = {'0.04', '0.09', '0.12', '0.16'};
dataset_labels = {'Static Input', 'Quenched Input', 'OU input'};
jseweights = {'0', '14','20','27'};
%%%%%%%%%%

%%%%%   %%%%%    %%%%%   %%%%%
cDir = [pwd '/'];

%%%%%   %%%%%    %%%%%   %%%%%
%%%%%   %%%%%    %%%%%   %%%%%

if external == 1
    %% PARSE ARGUMENTS
    %%%

    ddata = {dataOUexStat_27, dataH};
    inputtype  = {'OU', 'Quenched'};
    for cases = 1:numel(ddata)
        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        % Data cases
        casesexternal = ddata{cases};
        data1 = ddata{cases};
        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        P = parsePairs(varargin);
        checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

        % SETUP BASICS
        setPlotOpt('custom','path',cDir,'width',16,'height',9);
        Sep = '/';

        % PREPARE FIGURE
        figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
        DC = axesDivide(4,2,[0.09 0.09 0.75 0.75], .5, .45)';
        Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
        for i = 1:numel(DC)
            AH(i) = axes('Pos',DC{i}); hold on;
            set(0,'DefaultAxesTitleFontWeight','normal');
        end
        if P.Recompute
            LF_generateData(fname);
        end
        HF_setFigProps;

        limitval =1.0;
        Iapp = -1*limitval:0.1:limitval;
        idx = abs(Iapp)<=limitval;

        Sz = linspace(0.25,40,length(Iapp));%0.5*(1:21);
        corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
            0.8500, 0.3250, 0.0980; %p-p
            0.4660 0.6740 0.1880; %s-s %%som is green here
            0.4940, 0.1840, 0.5560]; %v-v
        colororder=corecolors+0.3;
        nullcolr = [0.6 0.6 0.6];
        Pops = {'E', 'PV', 'SOM', 'VIP'};
        Markers = {'o','o','o','o'};
        Linestyle = {'-','-','-'};

        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        % START PLOTTING
        for pp = 1:4
            pop = pp;
            axes(AH(pp));
            for k = 1:size(casesexternal,2)
                colr = colororder(pop,:)-0.15*(k-1); colr(colr>1)=1; colr(colr<0)=0;
                scatter(Iapp(idx),data1(k).CohPeakData.nuSim(idx,pop),Sz(idx),colr,'filled',Markers{1});
                plot(Iapp(idx),data1(k).CohPeakData.nuSim(idx,pop),'color',colr,'LineStyle',Linestyle{1},'linewidth',1.2);
            end
            plot(Iapp(idx),dataStat(1).CohPeakData.nuSim(idx,pop),'color',nullcolr,'linewidth',1.2)
            title(Pops{pop});
            ylabel([Pops{pop} ' Firing Rate (Hz)']);
            xlabel('input mean');
        end
        for pp = 1:4
            pop = pp;
            axes(AH(4+pp));
            for k = 1:size(casesexternal,2)
                colr = colororder(pop,:)-0.15*(k-1); colr(colr>1)=1; colr(colr<0)=0;
                scatter(data1(k).CohPeakData.nuSim(idx,pop),data1(k).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1});
                plot(data1(k).CohPeakData.nuSim(idx,pop),data1(k).CohPeakData.firstPeak(1,idx),'color',colr,'LineStyle',Linestyle{1},'linewidth',1.2);
            end
            plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',nullcolr,'linewidth',1.2);

            ylabel('E Max Coherence');
            xlabel([Pops{pop} ' Firing Rate (Hz)']);
            ylim([0 0.8]);
        end
        text(18,0.6,0,[sprintf('%s',inputtype{cases}),' var']);
       % ll = legend('',var_label{1},'',var_label{2});
        ll = legend('', var_label{1},'', var_label{2},'', var_label{3},'',var_label{4});
        ll.Box = 'off';
        ll.Position(1)= .8;
        ll.Position(2)= 0.18;
        ll.ItemTokenSize = [7,10];
        sst = sgtitle(sprintf('%s External Input to SOM',inputtype{cases}));
        sst.FontSize = 12;

        ffigname = [figname sprintf('%s',num2str(15+cases))];

        % SAVE FIGURE
        set(gcf, 'Renderer', 'painters');
        %savefig([cDir ffigname '.fig']);
        HF_viewsave('path',cDir,'name',ffigname,'view',P.View,'save',P.Save,'format','pdf','res',600);

        close all;
        clear ffigname;
    end
end

%%%%%%
if recurrent == 1
    %% PARSE ARGUMENTS
    dataStatnew(1) = dataStat(2); 
    dataStatnew(2) = dataStat(3); 
    dataStatnew(3) = dataStat(1); 

    ddata = {dataStatnew, dataOUrecE};
    inputtype  = {'Static Only', 'OU like E'};
    inputtypefile  = {'staticJse', 'OUlikeEJse'};

    jseWset = {{'14','22','27'}, jseweights};

    for cases = numel(ddata)
        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        % Data cases
        casesexternal = ddata{cases};
        data1 = ddata{cases};
        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        P = parsePairs(varargin);
        checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

        % SETUP BASICS
        setPlotOpt('custom','path',cDir,'width',16,'height',9);

        % PREPARE FIGURE
        figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
        DC = axesDivide(4,2,[0.09 0.09 0.75 0.75], .5, .45)';
        Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
        for i = 1:numel(DC)
            AH(i) = axes('Pos',DC{i}); hold on;
            set(0,'DefaultAxesTitleFontWeight','normal');
        end
        if P.Recompute
            LF_generateData(fname);
        end
        HF_setFigProps;

        limitval =1.0;
        Iapp = -1*limitval:0.1:limitval;
        idx = abs(Iapp)<=limitval;
        Sz = linspace(0.25,40,length(Iapp));%0.5*(1:21);
        corecolors = [0, 0.4470, 0.7410; %e-e
            0.8500, 0.3250, 0.0980; %p-p
            0.4660 0.6740 0.1880; %s-s %%som is green here
            0.4940, 0.1840, 0.5560]; %v-v
        colororder=corecolors+0.3;
        nullcolr = [0.6 0.6 0.6];

        Pops = {'E', 'PV', 'SOM', 'VIP'};
        Markers = {'o','o','o','o'};
        Linestyle = {'-','-','-'};

        %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%%%%
        % START PLOTTING
        jselabel = jseWset{cases};
        for pp = 1:4
            pop = pp;%subpop(pp)
            axes(AH(pp));
            for k = 1:size(casesexternal,2)
                colr = colororder(pop,:)-0.15*(k-1); colr(colr>1)=1; colr(colr<0)=0;
                scatter(Iapp(idx),data1(k).CohPeakData.nuSim(idx,pop),Sz(idx),colr,'filled',Markers{1})
                plot(Iapp(idx),data1(k).CohPeakData.nuSim(idx,pop),'color',colr,'LineStyle',Linestyle{1},'linewidth',1.2)
            end
            if cases~=1
                plot(Iapp(idx),dataStat(1).CohPeakData.nuSim(idx,pop),'color',nullcolr,'linewidth',1.2)
            end
            title(Pops{pop});
            ylabel([Pops{pop} ' Firing Rate (Hz)']);
            xlabel('static input');
        end

        for pp = 1:4%numel(subpop)
            pop = pp;
            axes(AH(4+pp));
            for k = 1:size(casesexternal,2)
                colr = colororder(pop,:)-0.15*(k-1); colr(colr>1)=1; colr(colr<0)=0;
                scatter(data1(k).CohPeakData.nuSim(idx,pop),data1(k).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
                plot(data1(k).CohPeakData.nuSim(idx,pop),data1(k).CohPeakData.firstPeak(1,idx),'color',colr,'LineStyle',Linestyle{1},'linewidth',1.2)
                text(0.9,0.78-0.12*k,jselabel{k},'unit','n','color',colr);
            end
            if cases~=1
                plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',nullcolr,'linewidth',1.2)
            end
            ylabel('E Max Coherence');
            xlabel([Pops{pop} ' Firing Rate (Hz)']);
            ylim([0 0.8]);
        end
        text(0.8,0.83,'J_{E\rightarrowSOM}','unit','n','color','k')   
        sst = sgtitle(sprintf('%s: Recurrent Input to SOM',inputtype{cases}));
        sst.FontSize = 12;

        ffigname = [figname '12'];

        % SAVE FIGURE
        set(gcf, 'Renderer', 'painters')
        %savefig([cDir ffigname '.fig'])
        HF_viewsave('path',cDir,'name',ffigname,'view',P.View,'save',P.Save,'format','pdf','res',600);

        close all;
    end
end

end