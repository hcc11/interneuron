function srr = poisson_distribution(rate,Ncell,Npop,T,N)
% srr generates spike train, based on genXspk code
% loops through each population of neurons generating spks at specified rates

Nsum= [0, cumsum(Ncell)];
srr=[];
Nsource=1;
for ii=1:Npop % each population is given a unique rate where spk trains are generated within that pop
    for ns=1:Nsource
        tempspikes=cumsum(-log(rand(1,round(rate(ii)*(Ncell(ii)/Nsource)*T*1.1)))/(rate(ii)*(Ncell(ii)/Nsource))); %generate spks at specificed rate
        tempspikes=tempspikes(tempspikes<T&tempspikes>0);  %filter spikes in proper time interval
        srr_temp=zeros(2,numel(tempspikes)); % generate temp matrix
        srr_temp(1,:)=tempspikes; % spktimes
        srr_temp(2,:)=ceil(rand(1,size(srr_temp,2))*Ncell(ii)/Nsource)+(ns-1)*Ncell(ii)/Nsource+Nsum(ii); %neuron ID %second term is 0
        srr=[srr,srr_temp]; % add to existing srr
    end
end

[~,J]=sort(srr(1,:)); % sort srr based on spk time
srr=srr(:,J);
srr=srr(:,srr(1,:)>0&srr(1,:)<=T); %% possible dupe from line 11 with tempspk filter

%%% plotting check - plots population spkrates match given rate input
% figure
% hold on;
% srr_rates = zeros(1,N);
% for nn=1:N
%     srr_rates(nn) = length(srr(1,srr(2,:) == nn))/(T);
% end
% xlabel('Neuron ID') 
% ylabel('Poisson Firing Rates') 
% scatter(1:N,srr_rates)
% hold off

end



%%%%%%% Old draft
%function [Iapp_full] = poisson_distribution(rate,current_mag, N,T,dt,Nt,Ncell,Npop)



%%% inputs
%Npop = ;%4; % number of populations (E, PV, SOM, VIP)
%N = 5e4; %total number of neurons in pop
%Ncell = 40;% number of neurons in pop %N_ratio/sum(N_ratio)* N;
%T=1000;
%dt=.05;  % time step (ms)
%Nt=T/dt;
%rIapp = 0.1;
%%%%
%Iapp_full = zeros(N,Nt);

%t = 0:dt:T-dt;

%Nsum= [0, cumsum(Ncell)];
%Nsource=1;%p.Nsource;

% tic;x
% sx=[];
% Nsource=1;
% for ii = 1:Npop
%     for ns=1:Nsource
%         tempspikes=cumsum(-log(rand(1,round(rate(ii)*(Ncell(ii)/Nsource)*T*1.1)))/(rate(ii)*(Ncell(ii)/Nsource))); %generate spk time
%         tempspikes=tempspikes(tempspikes<T&tempspikes>0); %filter spk time
%         sx_temp=zeros(2,numel(tempspikes));
%         sx_temp(1,:)=tempspikes;
%         sx_temp(2,:)=ceil(rand(1,size(sx_temp,2))*Ncell(ii)/Nsource)+(ns-1)*Ncell(ii)/Nsource+Nsum(ii); % generate neuron id
%         sx=[sx,sx_temp];
%     end
% end
%
% for ii = 1:Npop
%     for nn = Nsum(ii)+1:Nsum(ii+1)
%         spikeT = sx(1,sx(2,:) == nn); %grab the neuron I care about (nn) and the spike times
%         for tt = 2:length(t) % loop through all timesteps, at each dt bin count num of spks
%             numSpikes = length(spikeT(spikeT >= t(tt-1) & spikeT < t(tt))); %count all spikes in dt bin
%             Iapp_full(nn,tt-1) = current_mag(ii)*numSpikes; % at that dt bin, on nn neuron, numspks*current = iapp input on nn at that dt
%         end
%         %nn
%     end
% end
%
% %save('Iapp_spktrain')
%
% toc;