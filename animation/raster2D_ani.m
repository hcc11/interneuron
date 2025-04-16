function [varargout]=raster2D_ani(s0,t0,t1,Loc,Ncell,filename,data_folder,paramval)
%%%% initials movie to save
myVideo = VideoWriter(sprintf('%s%s_%0.2f.avi',data_folder,filename, paramval),'Motion JPEG AVI'); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
open(myVideo)

% Ne1: # of neurons per dimension 
    dta=2;
    %dta=2; % bin size for raster animation
    timea=t0:dta:t1; % Timeframes
%     timea=dta:dta:T; % Timeframes
%     spkduration=1; % Duration of each spike in raster
    fig1=figure('Visible','off');
    set(gcf,'color','w')
    set(gcf,'position',[350 100 550 450])
%     winsize = get(fig1,'Position');
%     winsize(1:2) = [0 0];
    numframes=numel(timea);
%     A=moviein(numframes,fig1,winsize);
    A(1:numframes) = struct('cdata', [],'colormap', []);
    for i=1:numframes
          Is=(s0(1,:)<=timea(i)& s0(1,:)>timea(i)-dta & s0(2,:)>0);
          idx = s0(2,Is); 
          Nsum = [0 cumsum(Ncell)]; 
          Npop = length(Ncell); 
          color = [0, 0.4470, 0.7410; %e-e  %%%% 
                0.8500, 0.3250, 0.0980; %p-p
                  0.4660 0.6740 0.1880; %s-s
                0.4940, 0.1840, 0.5560]; %v-v

          Size = [15 15 15 15]; 
          hold off 
          for pop = 1:Npop 
              tmp = idx(idx<=Nsum(pop+1) & idx>Nsum(pop)); 
              if isempty(tmp)
                  plot(-1,-1,'.','MarkerSize',Size(pop),'color',color(pop,:)) 
              else
                  plot(Loc(tmp,1), Loc(tmp,2),'.','MarkerSize',Size(pop),'color',color(pop,:))
              end
              hold on
          end
          hold off
          axis([0 1 0 1])
          xlabel('cell location X','fontsize',18)
          ylabel('cell location Y','fontsize',18)
          set(gca,'fontsize',15)
          set(gca,'xtick',[0 0.5 1])
          set(gca,'ytick',[0 0.5 1])
          
          title(sprintf('t=%d msec',round(timea(i))-t0))
         % legend('E','PV','SOM','VIP')
          pause(0.3)
          drawnow
          A(i)=getframe(gcf);
          writeVideo(myVideo,A(i));
    end
    if nargout==1 
        varargout{1}=A;
    end
    close(myVideo)