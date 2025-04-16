function SuppFig2(varargin)
%%%% Supplemental Figure 1

%%%%% Supp Fig 1 - shows the total average input (excitatory and
%%%%% inhibitory) with changes to static input

%%%%% Supp Fig 1 description - 2 rows x 4 cols
%%%% COLS - each input case; 
%%% col 1 - static input to E
%%% col 2 - static input to PV
%%% col 3 - static input to SOM
%%% col 4 - static input to VIP
%
%%%% ROWS 
%%% row 1 - (averaged) total input versus static input
%%% row 2 - (averaged) variance input versus static input


%%%% data loading

figname = 'SuppFig2';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',3); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',17,'height',17);
  
fnames = {'current_rate_Estatic.mat','current_rate_PVstatic.mat',...
    'current_rate_SOMstatic.mat','current_rate_VIPstatic.mat'}; 

outpath=cDir;

load('rate_est.mat'); 
for type = 1:4
     data(type) = load(fnames{type},'nuSim','I_mn','I_var','params','res');
end 

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,5,[0.06 0.06 0.9 0.9], .45, .3)';

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
limitval =1;
Iapp = data(1).params;
idx = abs(Iapp)<=limitval;
XLim = [-1*limitval limitval];

Sz = linspace(0.25,40,length(Iapp));%0.5*(1:21);
corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colororder=corecolors;

Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};


    
for type = 1:4
    iA = type;
    axes(AH(iA));
    for pop = 1:4
        plot(data(type).params(idx), data(type).I_mn(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
    end
    title(['Input to ' Pops{type}]);
    ylabel('avg $I_{e+i}$ current',Interpreter='latex');
    xlim([-1*limitval limitval]);
    set(gca,'xtick',-1:0.5:1)

    iA = type+4;
    axes(AH(iA));
    for pop = 1:4
        plot(data(type).params(idx), data(type).I_var(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
    end
    
    ylabel('var $I_{e+i}$ current',Interpreter='latex');
    xlim([-1*limitval limitval]);
    set(gca,'xtick',-1:0.5:1)

    iA = type+8;
    axes(AH(iA));
    for pop = 1:4
        plot(data(type).params(idx), data(type).nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
        plot(data(type).params(idx), rsim_EIF2(idx,pop,type),'o--','color',colororder(pop,:),'linewidth',0.5,'markersize',3)
%         plot(data(type).params(idx), r_est_selfcons(idx,pop,type),'o--','color',colororder(pop,:),'linewidth',0.5,'markersize',4)

    end
    xlim([-1*limitval limitval]);
    ylabel('firing rate (Hz)')
%     xlabel('static input');
set(gca,'xtick',-1:0.5:1)

    iA = type+12;
    axes(AH(iA));
    for pop = 1:4
        plot(data(type).params(idx), data(type).nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
        plot(data(type).params(idx), rsim_EIF3_2(idx,pop,type),'o--','color',colororder(pop,:),'linewidth',0.5,'markersize',3)

    end
    xlim([-1*limitval limitval]);
    ylabel('firing rate (Hz)')
    set(gca,'xtick',-1:0.5:1)
%     xlabel('static input');

    iA = type+16;
    axes(AH(iA));
    for pop = 1:4
        plot(data(type).params(idx), data(type).nuSim(idx,pop),'color',colororder(pop,:),'linewidth',1.2)
%         plot(data(type).params(idx), rsim_EIF(idx,pop,type),'o--','color',colororder(pop,:),'linewidth',0.5,'markersize',4)
        plot(data(type).params(idx), r_est_selfcons(idx,pop,type),'o--','color',colororder(pop,:),'linewidth',0.5,'markersize',3)

    end
    xlim([-1*limitval limitval]);
    ylabel('firing rate (Hz)')
    xlabel('static input');
    set(gca,'xtick',-1:0.5:1)
end

HF_setFigProps;

% % SAVE FIGURES
% savefig([outpath figname '.fig']);
 set(gcf, 'Renderer', 'painters');
 HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end

