function Figure6(varargin)

%%%%%% FIGURE 6 - OU input to SOM
%%% Figure 6 decription - 2 panels are generated. First panel plots changes
%%% of J_E-> SOM with OU like I_E->SOM input targetting SOM. SOM also has
%%% static applied to it.
%%% Plot 1: E max coh versus E FR
%%% note that grey plot in panel 1 is default static network case
%%% 
%%% Second panel corresponds to inset of paper figure (green plot here)
%%% with changes to the varuience of the OU like E input to SOM when
%%% J_E->SOM = 0 and SOM is targetted with static input.

close all

%%%%%   %%%%%    %%%%%   %%%%%
figname = 'Figure6';

%%%%%%%%%%
%%% load datasets
%%% static input
dataStat(1) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');

%%% OU replace E + static input to SOM: recurrent/Jse mods
dataOUrecE(1) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','0'),'CohPeakData');
dataOUrecE(2) = load(sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','14'),'CohPeakData');
dataOUrecE(3) = load('AllCohPwr_Avgs_resamp_OUtoSOMvar12_Jse20','CohPeakData');
dataOUrecE(4) = load( sprintf('AllCohPwr_Avgs_resamp_SOM_Jse%s_OUrec','27'),'CohPeakData');

%%% OU replace E + static input to SOM: recurrent/Jse mods
dataOUrecE_0(1) =  dataOUrecE(1);

datarecordSOM(1) = load('RecordSOM');
datareplaySOM(1) = load('ReplaySOM');

cDir = [pwd '/'];

%%%% %% default for comparison
%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);
setPlotOpt('custom','path',cDir,'width',19,'height',5);

% PREPARE FIGURE
fig = figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
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
tempcolor = corecolors(3,:)+0.15*[1,1,1];
colrshift = .2; pop = 1;

AH(1) = subplot(2, 8, [1 2 9 11]);
hold on
colrss = {tempcolor, [0.6 0.6 0.6], [0.2 0.2 0.2]};
style = {'-',':',':'};
var_label = {'0.0723', '0.04', '0.09'};
data = dataOUrecE_0; ss = 1;
colr = colrss{ss};
scatter(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
plot(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.5,'LineStyle',style{ss})
ylabel('E max coherence');
xlabel(sprintf('%s Firing Rate (Hz)',Pops{pop}));
ylim([0, 0.022]);
box off
AH(1).Position(1) = 0.07;
AH(1).Position(2) = 0.2;
AH(1).Position(3) = 0.15;
AH(1).Position(4) = 0.62;
hold off


%%%%%%%%%
%%% small quad figures 

%%% left column
AH(2) = subplot(2, 8, 3);
data = datareplaySOM;
colr = corecolors(1,:);
tmp = zeros(4,length(data.Iapp_range));
for tt = 1:data.Ntrial
    for pop = 1:4
        tmp(pop,:) = tmp(pop,:)+data.rate(pop,:,tt);
    end
end 
tmp = tmp/data.Ntrial;
plot(data(1).Iapp_range,tmp(1,:),'color',colr,'linewidth',1.5)
ylabel('rate');
box off
AH(2).Position(1) = 0.29;
AH(2).Position(2) = 0.62;
AH(2).Position(3) = 0.07;
AH(2).Position(4) = 0.25;

etmprate = tmp(1,:);

AH(4) = subplot(2, 8, 11);
colr = corecolors(1,:);
plot(data(1).Iapp_range,data(1).maxcohr(:,1),'color',colr,'linewidth',1.5)
ylabel('max Cohr');
xlabel('input');
box off
AH(4).Position(1) = AH(2).Position(1);
AH(4).Position(2) = AH(1).Position(2);
AH(4).Position(3) = AH(2).Position(3);
AH(4).Position(4) = AH(2).Position(4);

%%% right column
AH(3) = subplot(2, 8, 4);
data = datarecordSOM;
colr = corecolors(3,:);
plot(data(1).Iapp_range,mean(data(1).rate(3,:,:),3),'color',colr,'linewidth',1.5)
ylabel('rate');
box off
AH(3).Position(1) = 0.42;
AH(3).Position(2) = AH(2).Position(2);
AH(3).Position(3) = AH(2).Position(3);
AH(3).Position(4) = AH(2).Position(4);

AH(5) = subplot(2, 8, 12);
colr = corecolors(3,:);
plot(data(1).Iapp_range,data(1).maxcohr(:,3),'color',colr,'linewidth',1.5)
ylabel('max Cohr');
xlabel('input');
box off
AH(5).Position(1) = AH(3).Position(1);
AH(5).Position(2) = AH(4).Position(2);
AH(5).Position(3) = AH(2).Position(3);
AH(5).Position(4) = AH(2).Position(4);

%%%%%%
pop = 1;

AH(6) = subplot(2, 8, [5 6 13 14]);
%axes(AH(4)); %%%
hold on
colr = [0.6 0.6 0.6]; %%% default case
plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.7')
data = datareplaySOM; ss = 1;


replayIapp = data.Iapp_range;
idx1  = find(round(Iapp,2) == round(replayIapp(1),2));
idx2  = find(round(Iapp,2) == round(replayIapp(end),2));

newSz = linspace(Sz(idx1),Sz(idx2),length(replayIapp));


    colr = tempcolor-(ss-1)*[1 1 1]*colrshift;
    colr(colr>1) = 1;
    colr(colr<0) = 0;

    scatter(mean(data(1).rate(1,:,:),3),data(1).maxcohr(:,1),newSz,colr,'filled',Markers{1})
    plot(mean(data(1).rate(1,:,:),3),data(1).maxcohr(:,1),'color',colr,'linewidth',1.5)

ylabel('E max coherence');
xlabel(sprintf('%s Firing Rate (Hz)',Pops{pop}));
%ylim([0, 0.8]);
AH(6).Position(1) = 0.56;
AH(6).Position(2) = AH(1).Position(2);
AH(6).Position(3) = AH(1).Position(3);
AH(6).Position(4) = AH(1).Position(4);
hold off

AH(7) = subplot(2, 8, [7 8 15 16]);
%axes(AH(4)); %%%
hold on
colr = [0.6 0.6 0.6]; %%% default case
plot(dataStat(1).CohPeakData.nuSim(idx,pop),dataStat(1).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.7')

Jse_label = {'0', '14', '20', '27'};
data = dataOUrecE;
for ss = 1:size(data,2)
    colr = tempcolor-(ss-1)*[1 1 1]*colrshift;
    colr(colr>1) = 1;
    colr(colr<0) = 0;
    scatter(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),Sz(idx),colr,'filled',Markers{1})
    plot(data(ss).CohPeakData.nuSim(idx,pop),data(ss).CohPeakData.firstPeak(1,idx),'color',colr,'linewidth',1.5')
    if ss == 4
        text(0.5,0.8,'J_{E\rightarrowSOM}','unit','n','color','k')
    end
    text(0.7,0.75-0.1*ss,Jse_label{ss},'unit','n','color',colr);
end
ylabel('E max coherence');
xlabel(sprintf('%s Firing Rate (Hz)',Pops{pop}));
ylim([0, 0.8]);
AH(7).Position(1) = 0.77;
AH(7).Position(2) = AH(1).Position(2);
AH(7).Position(3) = AH(1).Position(3);
AH(7).Position(4) = AH(1).Position(4);
hold off

HF_setFigProps;

% SAVE FIGURES
set(gcf, 'Renderer', 'painters');
%savefig([cDir figname '.fig']);
HF_viewsave('path',cDir,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

clear figname
% close all

end



