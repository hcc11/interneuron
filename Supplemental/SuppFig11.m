function SuppFig11(varargin)
close all

%%% Supplemental Figure 6 - Ramping cases for bistability conditions
%%%% Output:
%%% Row 1: PV input
%%%   =>> Default Network: Jes = -120, Jps = -60
%%%     Col 1: E FR vs E Coh - Static Input
%%%     Col 2: Ramp Up E Trace and Ramp Down E Trace - indicating common
%%%                            overlap of input
%%%     Col 3: E FR vs E Coh - Average FR Across raming input time bins (5
%%%                            sec interval with +/-0.05 input change)
%%%   =>> Modified Netowrk: Jes = -60, Jps = -240
%%%     Col 4: E FR vs E Coh - Static Input
%%%     Col 5: Ramp Up E Trace and Ramp Down E Trace - indicating common
%%%                            overlap of input
%%%     Col 6: E FR vs E Coh - Average FR Across raming input time bins (5
%%%                            sec interval with +/-0.05 input change)
%%%
%%% Row 2: SOM input
%%%        Same as above but for SOM input.
%%%

%%%% SOM input
%%%% default
staticSOM(1) = load('AllCohPwr_Avgs_resamp_SOM_static','CohPeakData');
staticSOMdown(1) = load('SOMInputrampDown_Jes120Jps60');
staticSOMup(1) = load('SOMInputrampUp_Jes120Jps60');

%%%% mod circuit
staticmodSOM(1) = load('AllCohPwr_Avgs_resamp_Jes60_Jps240_somstatic','CohPeakData');
modSOMdown = load('SOMInputrampDOWN_Jes60Jps240');
modSOMup = load('SOMInputrampUP_Jes60Jps240');

%%% PV input
staticPV(1) = load('AllCohPwr_Avgs_resamp_PVstatic','CohPeakData');
staticPVdown(1) = load('PVInputrampDown_Jes120Jps60');
staticPVup(1) = load('PVInputrampUp_Jes120Jps60');  

%%%% mod circuit
staticmodPV(1) = load('AllCohPwr_Avgs_resamp_Jes60_Jps240_pvstatic','CohPeakData');
modPVdown = load('PVInputrampDOWN_Jes60Jps240');
modPVup = load('PVInputrampUP_Jes60Jps240');


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];

% setPlotOpt('custom','path',cDir,'cols',1,'height',17);
setPlotOpt('custom','path',cDir,'width',27,'height',10); 

inpath=[cDir ''];
outpath=[cDir ''];
Sep = '/';

% PREPARE FIGURE
fig = figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
if P.Recompute
    LF_generateData(fname);
end
HF_setFigProps;

Iapp = -1:0.1:1;
idx = abs(Iapp)<=1;

Sz = 1.5*(1:21)+5;
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v;

Pops = {'E', 'PV', 'SOM', 'VIP'};
% Markers = {'o','s','d','^'};
Markers = {'o','o','o','o'};


fig.Visible = 'off';

upcolr = [0.7 0.7 0.7];
downcolr = [0.3 0.3 0.3];

AH(1) = subplot(4, 6, [1 7]);
pop = 1; %% for max coherence data
hold on
scatter(staticPV(1).CohPeakData.nuSim(idx,pop),staticPV(1).CohPeakData.firstPeak(pop,idx),Sz(idx),corecolors(2,:),'filled');
plot(staticPV(1).CohPeakData.nuSim(idx,pop),staticPV(1).CohPeakData.firstPeak(pop,idx),'color',corecolors(2,:),'linewidth',1.2);
scatter(staticPV(1).CohPeakData.nuSim(10,pop),staticPV(1).CohPeakData.firstPeak(pop,10),Sz(10),corecolors(2,:),'filled',MarkerEdgeColor = 'b',LineWidth=2.0);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);

hold off

AH(1).Position(1) = AH(1).Position(1)-0.065;
AH(1).Position(2) = AH(1).Position(2)+0.05;
AH(1).Position(3) = AH(1).Position(3)-0.025;
AH(1).Position(4) = AH(1).Position(4)-0.09;
tt = text(7,0.1,0,'E Firing Rate (Hz)',Rotation=90,FontSize=9);
tt = text(-0.5,0.3,0,'PV input',Rotation=90,FontSize=10,Color=corecolors(2,:),FontWeight='bold');


%AH(1).Position = [ 0.1300    0.5482    0.1023    0.3768];
%%% rates 
%%% change colors and also have split in up and down
AH(2) = subplot(4, 6, 2);
hold on
t1 = 65000; t2 = 80000; %%75000 ms == static = -0.1
plot(t1:t2,staticPVup(1).RampCalc.re(t1:t2,1),'Color',upcolr);
xticklabels([]);

AH(2).Position(1) = AH(2).Position(1)-0.075;
AH(2).Position(2)= AH(2).Position(2)+0.012;%-0.02;
AH(2).Position(3) = AH(2).Position(3)+0.08;
AH(2).Position(4) = AH(2).Position(4)-0.02;
set(gca,'box','off');
tt = text(63000,26,0,'$J_{SOM\rightarrow E} = -120$ and $J_{SOM\rightarrow PV} = -60$',FontSize=12,FontWeight='bold',Interpreter='latex');
% AH(2).Position = [0.2645    0.7844    0.0993    0.1406];

AH(3) = subplot(4, 6, 8);
hold on
t1 = 85000; t2 = 100000; %%95000 ms == static = -0.1
plot(t1:t2,staticPVdown(1).RampCalc.re(t1:t2,1),'Color',downcolr);
%ylabel('E Firing Rate (Hz)');
xlabel('time (ms)');
xticklabels([]);

AH(3).Position(1) = AH(2).Position(1);
AH(3).Position(2)= AH(2).Position(2)-0.2;
AH(3).Position(3) = AH(2).Position(3);
AH(3).Position(4) = AH(2).Position(4);
set(gca,'box','off');

%%% ramp averages
AH(4) = subplot(4, 6, [3 9]);
hold on
plot(staticPVdown(1).RampCalc.nuSim_calcFRramp(5,:),staticPVdown.RampCalc.nuSim_calcFRramp(1,:),'Color',downcolr,LineWidth=1.4);
plot(staticPVup(1).RampCalc.nuSim_calcFRramp(5,:),staticPVup.RampCalc.nuSim_calcFRramp(1,:),'Color',upcolr,LineWidth=1.1);
ylabel('E Firing Rate (Hz)');
xlabel('input');
xlim([-1 1]);

hold off

AH(4).Position(1) = AH(4).Position(1)+0.005;%0.01;
AH(4).Position(2) = AH(1).Position(2);
AH(4).Position(3) = AH(1).Position(3);
AH(4).Position(4) = AH(1).Position(4);

annotation(fig,'rectangle',...
    [0.433,0.73,0.005,0.07],'FaceColor','r','FaceAlpha',.3,'EdgeColor','b');
annotation(fig,'rectangle',...
    [0.31,0.575,0.005,0.35],'FaceColor','r','FaceAlpha',.3,'EdgeColor','b');

x = [0.48 0.48];
y = [0.83 0.83];
annotation('textarrow',x,y,'String',{'input = -0.1'},'Color',corecolors(2,:),'HeadWidth',0,'HeadLength',0,'FontSize',8);

%  AH(4).Position =  [0.3991    0.5482    0.1023    0.3768]

%%%% next case
AH(5) = subplot(4, 6, [4 10]);
hold on
pop = 1; %% for max coherence data
scatter(staticmodPV(1).CohPeakData.nuSim(idx,pop),staticmodPV(1).CohPeakData.firstPeak(pop,idx),Sz(idx),corecolors(2,:),'filled');
plot(staticmodPV(1).CohPeakData.nuSim(idx,pop),staticmodPV(1).CohPeakData.firstPeak(pop,idx),'color',corecolors(2,:),'linewidth',1.2);
scatter(staticmodPV(1).CohPeakData.nuSim(17,pop),staticmodPV(1).CohPeakData.firstPeak(pop,17),Sz(17),corecolors(2,:),'filled',MarkerEdgeColor = 'b',LineWidth=2.0);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
hold off

AH(5).Position(1) = AH(1).Position(1)+0.48;
AH(5).Position(2) = AH(1).Position(2);
AH(5).Position(3) = AH(1).Position(3);
AH(5).Position(4) = AH(1).Position(4);
tt = text(14.75,0.1,0,'E Firing Rate (Hz)',Rotation=90,FontSize=9);

% AH(5).Position = [0.5336    0.5482    0.1023    0.3768];

%%% ramp averages
AH(8) = subplot(4, 6, [6 12]);
hold on
plot(modPVdown(1).RampCalc.nuSim_calcFRramp(5,:),modPVdown.RampCalc.nuSim_calcFRramp(1,:),'Color',downcolr,'linewidth',1.4);
plot(modPVup.RampCalc.nuSim_calcFRramp(5,:),modPVup.RampCalc.nuSim_calcFRramp(1,:),'Color',upcolr,'linewidth',1.1);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
xlim([-1 1]);
hold off

AH(8).Position(1) =  AH(4).Position(1)+0.49;
AH(8).Position(2) = AH(1).Position(2);
AH(8).Position(3) = AH(1).Position(3);
AH(8).Position(4) = AH(1).Position(4);

%AH(8).Position(1) = AH(8).Position(1)+0.01;
%AH(8).Position = [0.8027    0.5482    0.1023    0.3768];

%%% rates
%%% change colors and also have split in up and down
AH(6) = subplot(4, 6, 5);
t1 = 140000; t2 = 170000; %% 145000 ms == static = 0.6
plot(t1:t2,modPVup(1).RampCalc.re(t1:t2,1),'Color',upcolr);
%ylabel('E Firing Rate (Hz)');
%xlabel('up time');
xticklabels([]);

AH(6).Position(1) = AH(2).Position(1)+0.48;
AH(6).Position(2)= AH(2).Position(2);
AH(6).Position(3) = AH(2).Position(3);
AH(6).Position(4) = AH(2).Position(4);
set(gca,'box','off');
tt = text(135000,90,0,'$J_{SOM\rightarrow E} = -60$ and $J_{SOM\rightarrow PV} = -240$',FontSize=12,FontWeight='bold',Interpreter='latex');


AH(7) = subplot(4, 6, 11);
t1 = 40000; t2 = 70000; %% 45000 ms == static = 0.6
plot(t1:t2,modPVdown(1).RampCalc.re(t1:t2,1),'Color',downcolr);
%ylabel('E Firing Rate (Hz)');
xlabel('time (ms)');
xticklabels([]);

AH(7).Position(1) = AH(3).Position(1)+0.48;
AH(7).Position(2)= AH(3).Position(2);
AH(7).Position(3) = AH(2).Position(3);
AH(7).Position(4) = AH(2).Position(4);
set(gca,'box','off');

annotation(fig,'rectangle',...
    [0.955,0.68,0.005,0.15],'FaceColor','r','FaceAlpha',.3,'EdgeColor','b');
annotation(fig,'rectangle',...
    [0.7,0.575,0.005,0.35],'FaceColor','r','FaceAlpha',.3,'EdgeColor','b');

x = [0.985 0.985];
y = [0.86 0.86];
annotation('textarrow',x,y,'String',{'input = 0.6'},'Color',corecolors(2,:),'HeadWidth',0,'HeadLength',0,'FontSize',8);



%%%%%%% row 2 === SOM input
AH(9) = subplot(4, 6, [13 19]);
hold on
pop = 1; %% for max coherence data
scatter(staticSOM(1).CohPeakData.nuSim(idx,pop),staticSOM(1).CohPeakData.firstPeak(pop,idx),Sz(idx),corecolors(3,:),'filled');
plot(staticSOM(1).CohPeakData.nuSim(idx,pop),staticSOM(1).CohPeakData.firstPeak(pop,idx),'color',corecolors(3,:),'linewidth',1.2);
scatter(staticSOM(1).CohPeakData.nuSim(10,pop),staticSOM(1).CohPeakData.firstPeak(pop,10),Sz(10),corecolors(3,:),'filled',MarkerEdgeColor = 'b',LineWidth=2.0);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
hold off

AH(9).Position(1) = AH(1).Position(1);
AH(9).Position(2) = AH(1).Position(2)-0.42;
AH(9).Position(3) = AH(1).Position(3);
AH(9).Position(4) = AH(1).Position(4);
tt = text(9,0.1,0,'E Firing Rate (Hz)',Rotation=90,FontSize=9);
tt = text(-4.5,0.3,0,'SOM input',Rotation=90,FontSize=10,Color=corecolors(3,:),FontWeight='bold');

%AH(9).Position(1) =  AH(9).Position(1)-0.05;
% AH(9).Position = [ 0.1300    0.1100    0.1023    0.3768];

%%% ramp averages
%%%% Chanage colors!
AH(12) = subplot(4, 6, [15 21]);
hold on
plot(staticSOMdown(1).RampCalc.nuSim_calcFRramp(5,:),staticSOMdown.RampCalc.nuSim_calcFRramp(1,:),'Color',downcolr,LineWidth=1.4);
plot(staticSOMup(1).RampCalc.nuSim_calcFRramp(5,:),staticSOMup.RampCalc.nuSim_calcFRramp(1,:),'Color',upcolr,LineWidth=1.1);
ylabel('E Firing Rate (Hz)');
xlabel('input');
xlim([-1 1]);
hold off

AH(12).Position(1) = AH(4).Position(1);
AH(12).Position(2) = AH(4).Position(2)-0.42;
AH(12).Position(3) = AH(1).Position(3);
AH(12).Position(4) = AH(1).Position(4);
%AH(12).Position(1) = AH(12).Position(1)-0.05;
% AH(12).Position = [0.3991    0.1100    0.1023    0.3768];

%%% rates 
% %%% change colors and also have split in up and down
AH(10) = subplot(4, 6, 14);
t1 = 50000; t2 = 80000; %% 65000 ms == static = -0.1
plot(t1:t2,staticSOMup(1).RampCalc.re(t1:t2,1),'Color',upcolr);
%ylabel('E Firing Rate (Hz)');
%xlabel('up time');
xticklabels([]);

%AH(10).Position(1) = AH(10).Position(1)-0.05;
AH(10).Position(1) = AH(2).Position(1);
AH(10).Position(2)= AH(2).Position(2)-0.42;
AH(10).Position(3) = AH(2).Position(3);
AH(10).Position(4) = AH(2).Position(4);
set(gca,'box','off');


% AH(10).Position = [0.2645    0.3291    0.1023    0.1577];

AH(11) = subplot(4, 6, 20);
t1 = 80000; t2 = 110000; %% 95000 ms == static = -0.1
plot(t1:t2,staticSOMdown(1).RampCalc.re(t1:t2,1),'Color',downcolr);
%ylabel('E Firing Rate (Hz)');
xlabel('time (ms)');
xticklabels([]);

%AH(11).Position(1) = AH(11).Position(1)-0.05;
AH(11).Position(1) = AH(2).Position(1);
AH(11).Position(2) = AH(3).Position(2)-0.42;
AH(11).Position(3) = AH(2).Position(3);
AH(11).Position(4) = AH(2).Position(4);
set(gca,'box','off');

annotation(fig,'rectangle',...
    [0.433,0.355,0.005,0.09],'FaceColor',corecolors(3,:),'FaceAlpha',.3,'EdgeColor','b');
annotation(fig,'rectangle',...
    [0.28,0.16,0.005,0.35],'FaceColor',corecolors(3,:),'FaceAlpha',.3,'EdgeColor','b');
x = [0.497 0.497];
y = [0.433 0.433];
annotation('textarrow',x,y,'String',{'input = -0.1'},'Color',corecolors(3,:),'HeadWidth',0,'HeadLength',0,'FontSize',8);

AH(13) = subplot(4, 6, [16 22]);
hold on
pop = 1; %% for max coherence data
scatter(staticmodSOM(1).CohPeakData.nuSim(idx,pop),staticmodSOM(1).CohPeakData.firstPeak(pop,idx),Sz(idx),corecolors(3,:),'filled');
plot(staticmodSOM(1).CohPeakData.nuSim(idx,pop),staticmodSOM(1).CohPeakData.firstPeak(pop,idx),'color',corecolors(3,:),'linewidth',1.2);
scatter(staticmodSOM(1).CohPeakData.nuSim(6,pop),staticmodSOM(1).CohPeakData.firstPeak(pop,6),Sz(6),corecolors(3,:),'filled',MarkerEdgeColor = 'b',LineWidth=2.0);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
hold off

AH(13).Position(1) = AH(5).Position(1);
AH(13).Position(2) = AH(5).Position(2)-0.42;
AH(13).Position(3) = AH(1).Position(3);
AH(13).Position(4) = AH(1).Position(4);
tt = text(12,0.1,0,'E Firing Rate (Hz)',Rotation=90,FontSize=9);

%%% ramp averages
AH(16) = subplot(4, 6, [18 24]);
hold on
plot(modSOMdown(1).RampCalc.nuSim_calcFRramp(5,:),modSOMdown.RampCalc.nuSim_calcFRramp(1,:),'Color',downcolr,'linewidth',1.4);
plot(modSOMup.RampCalc.nuSim_calcFRramp(5,:),modSOMup.RampCalc.nuSim_calcFRramp(1,:),'Color',upcolr,'linewidth',1.1);
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
xlim([-1, 1]);
hold off

AH(16).Position(1) = AH(8).Position(1);
AH(16).Position(2) = AH(8).Position(2)-0.42;
AH(16).Position(3) = AH(1).Position(3);
AH(16).Position(4) = AH(1).Position(4);
% AH(16).Position = [0.8027    0.1100    0.1023    0.3768];

%%% rates
%%% change colors and also have split in up and down
AH(14) = subplot(4, 6, 17);
t1 = 40000; t2 = 90000; %% 55000 ms == static = -0.5
plot(t1:t2,modSOMup(1).RampCalc.re(t1:t2,1),'Color',upcolr);
%ylabel('E Firing Rate (Hz)');
%xlabel('up time');
xticklabels([]);

AH(14).Position(1) = AH(6).Position(1);
AH(14).Position(2) = AH(6).Position(2)-0.42;
AH(14).Position(3) = AH(2).Position(3);
AH(14).Position(4) = AH(2).Position(4);
set(gca,'box','off');

% AH(14).Position = [0.6682    0.3291    0.1023    0.1577];

AH(15) = subplot(4, 6, 23);
t1 = 140000; t2 = 190000; %% 155000 ms == static = -0.5
plot(t1:t2,modSOMdown(1).RampCalc.re(t1:t2,1),'Color',downcolr);
%ylabel('E Firing Rate (Hz)');
xlabel('time (ms)');
xticklabels([]);

AH(15).Position(1) = AH(7).Position(1);
AH(15).Position(2) = AH(7).Position(2)-0.42;
AH(15).Position(3) = AH(2).Position(3);
AH(15).Position(4) = AH(2).Position(4);
set(gca,'box','off');

annotation(fig,'rectangle',...
    [0.91,0.34,0.005,0.15],'FaceColor',corecolors(3,:),'FaceAlpha',.3,'EdgeColor','b');
annotation(fig,'rectangle',...
    [0.72,0.16,0.005,0.35],'FaceColor',corecolors(3,:),'FaceAlpha',.3,'EdgeColor','b');

x = [0.975 0.975];
y = [0.475 0.475];
annotation('textarrow',x,y,'String',{'input = -0.5'},'Color',corecolors(3,:),'HeadWidth',0,'HeadLength',0,'FontSize',8);

% Create rectangle
colr = ones(1,3)*0.85;
annotation(fig,'rectangle',...
    [0.0,0.94,1,0],LineStyle=':',LineWidth=1.5,color=colr);

annotation(fig,'rectangle',...
    [0.0,0.51,1,0],LineStyle=':',LineWidth=1.5,color=colr);

annotation(fig,'rectangle',...
    [0.0,0.09,1,0],LineStyle=':',LineWidth=1.5,color=colr);

annotation(fig,'rectangle',...
    [0.5,0,0,1],LineStyle=':',LineWidth=1.5,color=colr);

annotation(fig,'rectangle',...
    [0.03,0,0,1],LineStyle=':',LineWidth=1.5,color=colr);

annotation('textbox',...
    [0.2,0.0,0.1,0.1],...
    'String','Supercritical Hopf',...
    'FontSize',12,...
    'FontWeight','bold',...
    'Color',[0.7 0.7 0.7],...
    'EdgeColor','none');

annotation('textbox',...
    [0.7,0.0,0.1,0.1],...
    'String','Subcritical Hopf',...
    'FontSize',12,...
    'FontWeight','bold',...
    'Color',[0.7 0.7 0.7],...
    'EdgeColor','none');

x = [0.455 0.465];
y = [0.895 0.895];
annotation('textarrow',x,y,'String',{'(inc. input)','ramp up'},'Color',upcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);

x = [0.925 0.935];
y = [0.81 0.81];
annotation('textarrow',x,y,'String',{'ramp','up'},'Color',upcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);

x = [0.45 0.46];
y = [0.48 0.48];
annotation('textarrow',x,y,'String',{'ramp up'},'Color',upcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);

x = [0.922 0.932];
y = [0.315 0.315];
annotation('textarrow',x,y,'String',{'ramp', 'up'},'Color',upcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);



x = [0.467 0.467];
y = [0.71 0.71];
annotation('textarrow',x,y,'String',{'(dec.', 'input)'},'Color',downcolr,'HeadWidth',0,'HeadLength',0,'FontSize',8);
x = [0.445 0.435];
y = [0.64 0.64];
annotation('textarrow',x,y,'String',{'ramp', 'down'},'Color',downcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);


x = [0.94 0.93];
y = [0.64 0.64];
annotation('textarrow',x,y,'String',{'ramp','down'},'Color',downcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);

x = [0.445 0.435];
y = [0.22 0.22];
annotation('textarrow',x,y,'String',{'ramp','down'},'Color',downcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);

x = [0.94 0.93];
y = [0.26 0.26];
annotation('textarrow',x,y,'String',{'ramp','down'},'Color',downcolr,'HeadWidth',5,'HeadLength',2,'FontSize',8);




HF_setFigProps;

% SAVE FIGURES
% set(gcf, 'Renderer', 'opengl')
name = 'SuppFig11';
set(gcf, 'Renderer', 'painters')
%savefig([cDir name '.fig'])
HF_viewsave('path',cDir,'name',name,'view',P.View,'save',P.Save,'format','pdf','res',600);

close all;





end