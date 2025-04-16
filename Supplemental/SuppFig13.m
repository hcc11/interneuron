function SuppFig13(varargin) 


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);  

% SETUP BASICS
cDir =[pwd '/'];

setPlotOpt('custom','path',cDir,'width',18,'height',16); 

outpath=cDir;

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,3,[0.1 0.1 0.75 0.75], .35, .35)';% DC = axesDivide(4,1,[0.1 0.2 0.85 0.65], .5, 0.5)';

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
plottingpop = 2; % input to PV 
cases = {'default', 'sigmaSE', 'OtherSigmaSOM','AllsigmaSOM'}; 
figname = 'SuppFig13';
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
    tt = title(['Input to ' Pops{plottingpop} ': ' cases{k}],'FontSize',10);
    xlim(XLim);
%     ylim([0 maxyrates(k)+5]);
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
    pop=1; 
    for iii = 1:length(Iapp)
        if idx(iii) ==1
%             xx = (data(k).TrialAvgData.Cd_data_trialavg{iii}{pop,pop}...
%                 *data(k).TrialAvgData.Npair_data_trialavg{iii}{pop,pop})...
%                 /sum(data(k).TrialAvgData.Npair_data_trialavg{iii}{pop,pop});
            xx = mean(data(k).TrialAvgData.Cd_data_trialavg{1,iii}{1,1},2);
            tmp(:,count) = xx(1:51);
            if k==2
                 labels{1,count} = num2str(Iapp(iii));
            end
            count = count +1;
        end  
    end
    %%% Row 3: coherence versus frequency
    iA=k+8;
    axes(AH(iA));
    colors = cool(nnz(idx));
    freq = 0:50;
    hold on
    for iii = 1:size(cohplotcases{1,k},2)
        jjj = cohplotcases{1,k}(iii);
        plot(freq, tmp(:,jjj), 'Color',colors(jjj,:), 'LineWidth',1.15);  
    end
   if k == 1
    ylabel('E coherence','FontSize',10);
    if plottingpop == 4
        ylim([0 0.1]);
    end
   end
   if k == 2
       ll = legend(labels(cohplotcases{1,4}));
       ll.Title.String = 'Static Input';
       ll.Title.FontSize = 7;
       ll.Title.FontWeight = 'normal';
       ll.FontSize = 8;
       ll.FontWeight = 'normal';
       ll.Box = 'off';
       ll.Position(1) = 0.85;
       ll.ItemTokenSize = [7,10];
   end
    xlabel('Frequency (Hz)','FontSize',10);
    xlim([0 50]);
    xticks([0 25 50]);
    set(gca, 'FontSize',10);
    ylim([0 1])
end
HF_setFigProps;

% SAVE FIGURES
%savefig([outpath figname '.fig']);
set(gcf, 'Renderer', 'painters');
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);


end

