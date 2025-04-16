function SuppFig4(varargin)


figname = 'SuppFig4';

%% PARSE ARGUMENTS
P = parsePairs(varargin);
checkField(P,'FIG',2); checkField(P,'Save',1);  checkField(P,'View',1);  checkField(P,'Recompute',0);

% SETUP BASICS
cDir = [pwd '/'];
setPlotOpt('custom','path',cDir,'width',15,'height',10);
  
outpath=cDir;
load('respgain_deltainput.mat'); 


% PREPARE FIGURE
figure(P.FIG); clf; set(P.FIG,FigOpt{:}); HF_matchAspectRatio;
DC = axesDivide(3,2,[0.1 0.1 0.8 0.75], .45, .3)';

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


Pops = {'E', 'PV', 'SOM', 'VIP'};
Markers = {'o','o','o','o'};
Ylims = {[0 12],[0 2]}; 

Tseg = 100; 
t=(-Tseg/2+step/2):step:Tseg/2;
pid_range = 6:21; 
pid_plotrange = 6:3:21; 
colors = cool(length(pid_plotrange));  
Iapp = -1:0.1:1;  
for pop = 1 : 2
    iA = 1+(pop-1)*3;
    axes(AH(iA));
    for pp =1:length(pid_plotrange)
        pid = pid_plotrange(pp);
        y = res_deltastim(pid,pop).rate_ave- mean(res_deltastim(pid,pop).rate_ave(1:99));
        plot(t, y,'color',colors(pp,:),'linewidth',.5);
        text(.8,1-pp*0.1,sprintf('%.1f',Iapp(pid)),'unit','n','color',colors(pp,:))
    end
    text(.65,1,'input to PV','unit','n','color','k')
    xlabel('time (ms)')
    ylabel('rate (Hz)')

    a = zeros(length(pid_plotrange),1);
    tau = zeros(length(pid_plotrange),1);
    for pp =1:length(pid_range)
        pid = pid_range(pp);
        a(pp) = res_deltastim(pid,pop).f1.a;
        tau(pp) = res_deltastim(pid,pop).f1.tau;
    end
    iA = 2+(pop-1)*3;
    axes(AH(iA));
    plot(Iapp(pid_range),a,'o-','linewidth',1,'markersize',5)
    ylabel('a')
    xlabel('static input')

    iA = 3+(pop-1)*3;
    axes(AH(iA));
    plot(Iapp(pid_range),tau,'o-','linewidth',1)
    ylabel('\tau')
    ylim(Ylims{pop})
    xlabel('static input')
end


HF_setFigProps;

% % SAVE FIGURES
% savefig([outpath figname '.fig']);
 set(gcf, 'Renderer', 'painters');
 HF_viewsave('path',outpath,'name',figname,'view',P.View,'save',P.Save,'format','pdf','res',600);

end

