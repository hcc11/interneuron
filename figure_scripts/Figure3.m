function Figure3(varargin) 

%%%% Figure 3 displays each case of static input applied to one population
% average population rates 

%%%%% Figure description: 3 rows x 4 cols
%%%% COLS: each input case (E, PV, SOM, or VIP - in this column order)
%
%%%% ROWS: averaged (across trials) data for each poulation across static input for the 4
%%%% input cases
%%% row 1: FR
%%% row 2: max coherence 
%%% row 3: powerspectrum of E population coherence as function of frequency,
%%% each line is a single value of static input


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir =[pwd '/'];

setPlotOpt('custom','path',cDir,'width',18,'height',14); 
outpath=[cDir];
Sep = '/';

figname = 'Figure3';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,3,[0.06 0.1 0.85 0.85], .5, .5)';% DC = axesDivide(4,1,[0.1 0.2 0.85 0.65], .5, 0.5)';

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
data(1) = load('AllCohPwr_Avgs_resamp_Estatic','CohPeakData');
data(2) = load('AllCohPwr_Avgs_resamp_PVstatic','CohPeakData');
data(3) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');
data(4) = load('AllCohPwr_Avgs_resamp_VIPstatic','CohPeakData');

data_cohr(1) = load('AllCohPwr_Avgs_resamp_Estatic','TrialAvgData'); 
data_cohr(2) = load('AllCohPwr_Avgs_resamp_PVstatic','TrialAvgData'); 
data_cohr(3) = load('AllCohPwr_Avgs_resamp_SOM_static','TrialAvgData'); 
data_cohr(4) = load('AllCohPwr_Avgs_resamp_VIPstatic','TrialAvgData'); 

limitvalue = 0.7;
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

maxyrates = [30 20 25 80];

%%% for powerspectrum
%%%% only plotting a few static input cases in row 3 to avoid overcrowding
%%%% each panel
freq = 0:500; 
cohplotcases = {[1 4 8 12 15],[1 4 8 12 15],[1 4 8 12 15],[1 4 8 12 15]};
caseshere = sort(unique([1:2:15, 8]));
cohplotcases = {caseshere, caseshere, caseshere, caseshere};

for k = 1:4
    %%% Row 1: FR versus static input
    iA=k;
    axes(AH(iA));
    for pop=1:4
        plot(Iapp(idx), data(k).CohPeakData.nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
        if k ==4
            text(1.07,1.0-0.15*pop,Pops{pop},'unit','n','color',colororder(pop,:),'FontSize',8);
        end
        if k == 1
            ylabel('Firing rate (Hz)','FontSize',10);
        end
        xlabel('Static Input','FontSize',10);
    end
    tt = title(['Input to ' Pops{k}],'FontSize',10);
    xlim(XLim);
    ylim([0 maxyrates(k)+5]);
    set(gca, 'FontSize',10);
end

for k = 1:4
    %%% Row 2: avg maximum coherence versus static input
    iA=k+4;
    axes(AH(iA));
    for pop=1:4
        plot(Iapp(idx), data(k).CohPeakData.firstPeak(pop,idx),'color',colororder(pop,:),'linewidth',1.2)
    end
   if k == 1
    ylabel('Max coherence','FontSize',10);
   end
    xlabel('Static Input','FontSize',10);
    xlim(XLim);
    ylim([0 1]);
    set(gca, 'FontSize',10);
end

labels = cell(1,nnz(idx));
for k = 1:4
    tmp = zeros(51,nnz(idx));
    count = 1;
    for iii = 1:length(Iapp)
        if idx(iii) ==1
            xx = mean(data_cohr(k).TrialAvgData.Cd_data_trialavg{1,iii}{1,1},2);
            tmp(:,count) = xx(1:51);
            if k==4
                 labels{1,count} = num2str(Iapp(iii));
            end
            count = count +1;
        end  
    end
    %%% Row 3: coherence versus frequency
    iA=k+8;
    axes(AH(iA));
    colors = cool(nnz(idx));
    freq = 1:51;
    hold on
    for iii = 1:size(cohplotcases{1,k},2)
        jjj = cohplotcases{1,k}(iii);
            plot(freq, tmp(:,jjj), 'Color',colors(jjj,:), 'LineWidth',1.15);
       if k < 4
           ylim([0 1]);
       elseif k == 4 
       end

    end
   if k == 1
    ylabel('E coherence','FontSize',10);
   end
   if k == 4
       ll = legend(labels(cohplotcases{1,4}));
       ll.Title.String = 'Static Input';
       ll.Title.FontSize = 7;
       ll.Title.FontWeight = 'normal';
       ll.FontSize = 8;
       ll.FontWeight = 'normal';
       ll.Box = 'off';
       ll.Position(1) = 0.915;
       ll.ItemTokenSize = [7,10];
   end
    xlabel('Frequency (Hz)','FontSize',10);
    xlim([0 50]);
    xticks([0 25 50]);
    set(gca, 'FontSize',10);
end
HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig']);
set(gcf, 'Renderer', 'painters');
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);
close all

end

