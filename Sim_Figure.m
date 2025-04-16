function Sim_Figure(filename,T, input, popwihtinput, varargin)
filedata = load(filename);
addpath(genpath([pwd '/BE_figTools/']));


%% PARSE ARGUMENTS
% SETUP BASICS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

cDir =[pwd '/'];
setPlotOpt('custom','path',cDir,'width',6,'height',11);
outpath=[cDir];
Sep = '/';

figname = 'TestFig';


% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(1,[2 0.6 0.8 2] ,[0.2 0.1 0.65 0.85], .05, [0.5 0.5 1 0.5])';%[0.05 0.5 0.55 0.55])';
%DC = axesDivide(3 ,2 ,[0.1 0.2 0.85 0.6], 0.5, 0.5)';


Labels = {'A','B','C','D','E','F'}; LdPos = [-0.1,0.07];
for i = 1:numel(DC)
    AH(i) = axes('Pos',DC{i}); hold on;
    set(0,'DefaultAxesTitleFontWeight','normal');
end
if P.Recompute
    LF_generateData(fname);
end

% START PLOTTING
%%% plotting esthetics
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors;
Pops = {'E', 'PV', 'SOM', 'VIP'};
corrd_labels = {'E-E', 'PV-PV', 'SOM-SOM', 'VIP-VIP'};

freq = 0:500; %range of freq
daxis = 1/80*(1:2:40); %data_cohrIndivid.IndividData.CohPwrSpec_Calc{1,1}.daxis;
Tw = 8; % small window only filter over 850 seconds for plot
t2 = T-200; %snapshot of time simulation;
t1 = t2-1000;

N = 5e4;
Ncell = [40000 4000 4000 2000];
Nsum = [ 0 cumsum(Ncell)];
binsize = 1;
samplesize = 500;

% Row 1 : spk rasters
iA = 1;
axes(AH(iA));
spkdata = filedata.s1;
initialT = 500; finalT = T;
time = initialT:binsize:finalT;

popcountinit = 1;
for pop = 1:4
    nPOPcells = Ncell(pop);
    Popspks = spkdata(:,spkdata(2,:) <= Nsum(pop+1) & spkdata(2,:) > Nsum(pop));
    Popspks = Popspks(:, Popspks(1,:)>t1 & Popspks(1,:)<=t2);
    Popratio = (nPOPcells/50000)*samplesize;

    skips = 1;
    num_Popneurons_firing = size(unique(Popspks(2,:)),2); % number of neurons that are firing
    Popneurons_firing = unique(Popspks(2,:)); % Unique neurons IDs that are firing
    popcountfin = popcountinit + Popratio;

    if num_Popneurons_firing>= Popratio
        % rng(24);
        samplePopspks = randsample(Popneurons_firing,Popratio,false);
        hold on
        for ii = popcountinit:skips:popcountfin-1
            Popsamplespktimes = Popspks(1,Popspks(2,:) == samplePopspks(ii+1-popcountinit));
            scatter(Popsamplespktimes,ones(1,length(Popsamplespktimes))*ii,2,corecolors(pop,:),'filled');
        end

    elseif num_Popneurons_firing < Popratio && num_Popneurons_firing ~=0
        samplePopspks = randsample(Popneurons_firing,num_Popneurons_firing,false);
        hold on
        for ii = popcountinit:skips:length(samplePopspks)
            Popsamplespktimes = Popspks(1,Popspks(2,:) == samplePopspks(ii));
            scatter(Popsamplespktimes,ones(1,length(Popsamplespktimes))*ii,2,corecolors(pop,:),'filled');
        end
    end
    popcountinit = popcountfin;
    clear Popratio Popneurons_firing num_Popneurons_firing Popspks

end
% xlabel('Time (ms)','FontSize',11)
xlim([t1 t2])
ylim([0 samplesize])
% set(gca, 'XTick',[t1 ceil(t2+t1)/2 t2],'YTick',[],'FontSize',10)
set(gca, 'XTick',[],'YTick',[])
ylabel('Neuron','FontSize',11)
for pops = 1:4
    text(1.05,0.4+0.12*pops,[Pops{pops}],'unit','n','color',colororder(pops,:), 'FontSize',8);
end
title(sprintf('Static to %s, input = %0.2f',Pops{popwihtinput},input));


%%% rate traces
iA = 2;
axes(AH(iA));
re = filedata.corrCalc.reUnsmoothed;
re_smoothed=imfilter(re,ones(Tw,1)/Tw); %smooth rate

hold on
for pop=1:4
    plot((t1:t2), re_smoothed(t1:t2,pop)','color',colororder(pop,:),'linewidth',1.2)
end
text(-0.2,0.05,'Rate (Hz)','unit','n','Rotation',90,'FontSize',10);
xlim([t1 t2])
xlabel('Time (ms)','FontSize',11)
set(gca, 'XTick', [t1 ceil(t2+t1)/2 t2],'FontSize',10)
hold off
clear re

Tburn = 200;
re_each =hist(spkdata(2,spkdata(1,spkdata(2,:)>0)>Tburn),1:5e4)/(T-Tburn)*1e3;
mean_re = [mean(re_each(Nsum(1)+1:Nsum(2))),mean(re_each(Nsum(2)+1:Nsum(3))),...
    mean(re_each(Nsum(3)+1:Nsum(4))), mean(re_each(Nsum(4)+1:Nsum(5)))];
std_re = [std(re_each(Nsum(1)+1:Nsum(2))),std(re_each(Nsum(2)+1:Nsum(3))), ...
    std(re_each(Nsum(3)+1:Nsum(4))), std(re_each(Nsum(4)+1:Nsum(5)))];
stderror_ofmean_re = std_re./sqrt(Ncell);

%%% Row 3: bar plot for average rates
iA =3;
axes(AH(iA));
b= bar(filedata.nuSim,'FaceColor','flat');
b.CData = corecolors;
box off
hold on
for iii = 1:numel(b.XEndPoints)
    xtips = b.XEndPoints(iii);
    ytips = b.YEndPoints(iii);
    errorbar(xtips,ytips,stderror_ofmean_re(iii),'k.','LineWidth',1.2);
end
hold off

xticks([1:4]);
xticklabels({'E','PV','SOM','VIP'});
ax = gca;
ax.XAxis.FontSize = 11;
text(-0.2,0.0,'Rate (Hz)','unit','n','Rotation',90,'FontSize',10);
numsim = round(filedata.nuSim,2);
for pops = 1:4
    text(0.9,1.25-0.27*pops,[num2str(numsim(pops),2)],'unit','n','color',colororder(pops,:), 'FontSize',8)
    %        text(1.05,0.9-0.18*pops,['sem:' num2str(round(stderror_ofmean_re(pops),4))],'unit','n','color',colororder(pops,:), 'FontSize',10)
end

%%% Row 4: corr_d vs dist
iA =4;
axes(AH(iA));
hold on;
for pop = 1:4
            plot(filedata.corrCalc.daxis,filedata.corrCalc.C_d{pop,pop},'color',colororder(pop,:),'LineWidth',1.2);
            text(0.8,0.9-0.12*pop,[corrd_labels{pop}],'unit','n','color',colororder(pop,:), 'FontSize',8);
end
xlabel('Pairwise Distance (a.u.)');
ylim([-0.05 1]);
xlim([0 0.5]);
xticks([0:0.1:0.5]);
ylabel('Spike Count Correlation');




HF_setFigProps;

% SAVE FIGURES
set(gcf, 'Renderer', 'painters')
P.Save = 1;
%savefig([outpath figname '.fig'])
HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

% close all
end