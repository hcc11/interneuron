function Figure5(varargin)
close all

%%%% Figure 5 description: 4 rows x 4 cols
%%% each plot is the (y,x): E max coherence x E FR
%%% plots the cases of Pv static input and Som static input across modified
%%% circuit strengths:
%%% ROW 1: (neg) increasing values of Jes = {-60, -120, -240} with Jsp=0
%%%
%%% ROW 2: (neg) increasing values of Jps = {0, -120, -240, -420} woth Jes
%%% = -120
%%%
%%% COL 1: Static input to PV; COL 2: Static input to SOM


figname = 'Figure5';

%%%%%%%%%%%%%
% % START PLOTTING
% % default param Jes = -120, Jps = -60
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


%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',4); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',10,'height',9); 

Sep = '/';

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(2,2,[0.1 0.1 0.85 0.85], [0.75 0.0], 0.4)';

Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
     AH(i) = axes('Pos',DC{i}); hold on;
     DC{i};
    set(0,'DefaultAxesTitleFontWeight','normal');
end
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
corecolors(2,:) = corecolors(2,:)+0.2*[1,1,1];
corecolors(3,:) = corecolors(3,:)+0.2*[1,1,1];

Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};

for jj = 1:7
    if jj >1 && jj<5
        data_somJps0(jj-1) = load(datafnamesSOM{jj},'CohPeakData');
        data_pvJps0(jj-1) = load(datafnamesPV{jj},'CohPeakData');
    end
    if jj <5
        data_somJes60(jj) = load(datafnamesSOM{2,jj},'CohPeakData');
        data_pvJes60(jj) = load(datafnamesPV{2,jj},'CohPeakData');
    end
    if jj <7
        data_pvJes120(jj)=load(datafnamesPV{3,jj},'CohPeakData');
        data_somJes120(jj)=load(datafnamesSOM{3,jj},'CohPeakData');
    end
    data_pvJes240(jj)=load(datafnamesPV{4,jj},'CohPeakData');
    data_somJes240(jj)=load(datafnamesSOM{4,jj},'CohPeakData');
end

Jes_label = {'-60', '-120', '-240'};
Jps_label = {'  0','-60', '-120', '-240', '-360', '-420'};

iA = 1;
axes(AH(iA));
colrshift = 0.25; %%%%% changing J S->E/ J S->PV=0 (PV)
pop = 1; %% for max coherence data
for jj = size(data_pvJps0,2):-1:1
    tmp=data_pvJps0(jj).CohPeakData.nuSim(:,pop)<0.2;
    data_pvJps0(jj).CohPeakData.firstPeak(pop,tmp)=0;
 
    colr = corecolors(2,:)-(jj-1)*[1 1 1]*colrshift;
                colr(colr>1) = 1;
                colr(colr<0) = 0;
    scatter(data_pvJps0(jj).CohPeakData.nuSim(idx,pop),data_pvJps0(jj).CohPeakData.firstPeak(pop,idx),Sz(idx),colr,'filled');
    plot(data_pvJps0(jj).CohPeakData.nuSim(idx,pop),data_pvJps0(jj).CohPeakData.firstPeak(pop,idx),'color',colr,'linewidth',1.2);
    
    if jj == 1
            text(1.0,0.55,'$J_{SOM \rightarrow E}$','Interpreter','latex','unit','n','color','k','FontSize', 10);
    end
    fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
    text(1.0,0.65-0.12*(jj+1),Jes_label{jj},'unit','n','color',fontcolr,'FontSize', 8);
end
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
xlim([0 8]);
tt = title('$J_{SOM \rightarrow PV} = 0$','Interpreter','latex');
tt.Position(1) = 10;
tt.FontSize = 12;

iA = 3;
axes(AH(iA)); %%%%% changing J S->PV/ J S->E=-120 (PV)
pop = 1;%% for max coherence data
indexvals = [1, 3, 4, 6];
for  ff = 1:numel(indexvals)
    subdata_pvJes120(ff) = data_pvJes120(indexvals(ff)); 
    subdata_somJes120(ff) = data_somJes120(indexvals(ff));
end
subJps_label =  Jps_label(indexvals);

for jj = size(subdata_pvJes120,2)-1:-1:1
    tmp=subdata_pvJes120(jj).CohPeakData.nuSim(:,pop)<0.2;
    subdata_pvJes120(jj).CohPeakData.firstPeak(pop,tmp)=0;

    colr = corecolors(2,:)-(jj-1)*[1 1 1]*colrshift;
                colr(colr>1) = 1;
                colr(colr<0) = 0;

    scatter(subdata_pvJes120(jj).CohPeakData.nuSim(idx,pop),subdata_pvJes120(jj).CohPeakData.firstPeak(pop,idx),Sz(idx),colr,'filled');
    plot(subdata_pvJes120(jj).CohPeakData.nuSim(idx,pop),subdata_pvJes120(jj).CohPeakData.firstPeak(pop,idx),'color',colr,'linewidth',1.2);
    if jj == 1
            text(1.0,0.75,'$J_{SOM \rightarrow PV}$','Interpreter','latex','unit','n','color','k','FontSize', 10);
    end
    fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
    text(1.0,0.75-0.12*jj,subJps_label{jj},'unit','n','color',fontcolr,'FontSize',8);
end

    jj = size(subdata_pvJes120,2);
    tmp=subdata_pvJes120(jj).CohPeakData.nuSim(:,pop)<0.2;
    subdata_pvJes120(jj).CohPeakData.firstPeak(pop,tmp)=0;
    colr = corecolors(2,:)-(jj-1)*[1 1 1]*colrshift;
                colr(colr>1) = 1;
                colr(colr<0) = 0;
    scatter(subdata_pvJes120(jj).CohPeakData.nuSim(idx,pop),subdata_pvJes120(jj).CohPeakData.firstPeak(pop,idx),Sz(idx),colr,'filled');
    plot(subdata_pvJes120(jj).CohPeakData.nuSim(idx,pop),subdata_pvJes120(jj).CohPeakData.firstPeak(pop,idx),'color',colr,'linewidth',1.2);

    fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
    text(1.0,0.75-0.12*jj,subJps_label{jj},'unit','n','color',fontcolr,'FontSize',8);

xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
xlim([0 8]);
tt = title('$J_{SOM \rightarrow E} = -120$','Interpreter','latex');
tt.Position(1) = 10;
tt.FontSize = 12;

iA = 2;
axes(AH(iA));%%%%% changing J S->E/ J S->PV=0 (SOM)
pop = 1;%% for max coherence data
for jj = size(data_somJps0,2):-1:1

    tmp=data_somJps0(jj).CohPeakData.nuSim(:,pop)<0.2;
    data_somJps0(jj).CohPeakData.firstPeak(pop,tmp)=0;
    
    colr = corecolors(3,:)-(jj-1)*[1 1 1]*colrshift;
                colr(colr>1) = 1;
                colr(colr<0) = 0;

   % colr = cc(colrrr(jj),:);
    scatter(data_somJps0(jj).CohPeakData.nuSim(idx,pop),data_somJps0(jj).CohPeakData.firstPeak(pop,idx),Sz(idx),colr,'filled');
    plot(data_somJps0(jj).CohPeakData.nuSim(idx,pop),data_somJps0(jj).CohPeakData.firstPeak(pop,idx),'color',colr,'linewidth',1.2);
    fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
    text(-0.5,0.65-0.12*(jj+1),Jes_label{jj},'unit','n','color',fontcolr,'FontSize',8);
end
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
xlim([0 8]);

%%%%
iA = 4;
axes(AH(iA));%%%%% changing J S->PV/ J S->E=-120 (SOM)
pop = 1;
for jj = size(subdata_somJes120,2):-1:1
    tmp=subdata_somJes120(jj).CohPeakData.nuSim(:,pop)<0.2;
    subdata_somJes120(jj).CohPeakData.firstPeak(pop,tmp)=0;
    colr = corecolors(3,:)-(jj-1)*[1 1 1]*colrshift;
                colr(colr>1) = 1;
                colr(colr<0) = 0;
    scatter(subdata_somJes120(jj).CohPeakData.nuSim(idx,pop),subdata_somJes120(jj).CohPeakData.firstPeak(pop,idx),Sz(idx),colr,'filled');
    plot(subdata_somJes120(jj).CohPeakData.nuSim(idx,pop),subdata_somJes120(jj).CohPeakData.firstPeak(pop,idx),'color',colr,'linewidth',1.2);
    fontcolr = colr; %nullcolor+(kk-1).*[1 1 1]*colrshift;
    text(-0.5,0.75-0.12*jj,subJps_label{jj},'unit','n','color',fontcolr,'FontSize',8);
end
xlabel('E Firing Rate (Hz)');
ylabel('E max coherence');
ylim([0 1]);
xlim([0 8]);

HF_setFigProps;

% SAVE FIGURES
set(gcf, 'Renderer', 'painters')
%savefig([cDir figname '.fig'])
HF_viewsave('path',cDir,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);



end