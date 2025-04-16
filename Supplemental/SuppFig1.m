function SuppFig1(varargin)

%%%% data loading
figname = 'SuppFig1';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',17,'height',5);
  
fnames = {'current_rate_Estatic.mat','current_rate_PVstatic.mat',...
    'current_rate_SOMstatic.mat','current_rate_VIPstatic.mat'}; 

inpath=[cDir];
outpath=cDir;
load(fnames{1}); 

data(1) = load('FIcurve_Ie-0d2.mat');
data(2) = load('FIcurve_Ie0.mat');
data(3) = load('FIcurve_Ie0d4.mat');

% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(4,1,[0.05 0.15 0.92 0.75], .45, .3)';

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

corecolors = [0, 0.4470, 0.7410; %e-e  %%%% first four default matlab colros
    0.8500, 0.3250, 0.0980; %p-p
    0.4660 0.6740 0.1880; %s-s %%som is green here
    0.4940, 0.1840, 0.5560]; %v-v
colors=corecolors;

Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};
Xlims = {[0 1.5],[0 1.5],[-0.5 1.]}; 

    
for k = 1:3
    iA = k;
    axes(AH(iA));
    for pop =1:4
        cell_idx = (1:500)+(pop-1)*500;
%          scatter(res(data(k).pid).current_mn(cell_idx),res(data(k).pid).rate(cell_idx),2,colors(pop,:))
        plot(res(data(k).pid).current_mn(cell_idx),res(data(k).pid).rate(cell_idx),'.','color',colors(pop,:),'markersize',2)
        plot(data(k).FIcurve(pop).mu0, data(k).FIcurve(pop).rsim_vec, '-', 'LineWidth', 1,'color',colors(pop,:));
    if k==1 
        text(0.1,1-0.1*pop,Pops{pop},'unit','n','color',colors(pop,:))
    end
    end 
xlim(Xlims{k})
ylim([0 45]) 
    title(sprintf('Input to E: %.2g', params(data(k).pid)));
    ylabel('firing rate (Hz)');
    xlabel('mean current')

end

iA =4; 
axes(AH(iA));
load('FI_EIF.mat');
colors = cool(length(I_var_range)); 
for kk = 1: length(I_var_range)
    plot(FIcurve(kk).mu0, FIcurve(kk).rsim_vec, '-', 'LineWidth', 1,'color',colors(kk,:)); % Simulated firing rates
    text(0.1, 0.9-0.1*kk,sprintf('%.1f',I_var_range(kk)),'unit','n','color',colors(kk,:))
end 
text(0.1, 1,'I_{var}','unit','n','color','k')
ylabel('firing rate (Hz)');
xlabel('mean current')
xlim([-0.4 1.2])

HF_setFigProps;

% % SAVE FIGURES
% savefig([outpath figname '.fig']);
 set(gcf, 'Renderer', 'painters');
 HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end

