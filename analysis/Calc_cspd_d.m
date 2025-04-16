function Calc_cspd_d(data_folder, filename, samplesize, resampleCounts)

%%%%minActiveNum = samplesize*resampleCounts;

original_files = dir(filename);
file = load([data_folder original_files.name]);
matlabdata_filename = filename;

nametosave = extractBefore(filename,'.mat');

s1 = file.s1;
Npop = length(file.param.Ncell);
Nsum=[0 cumsum(file.param.Ncell)];
N = sum(file.param.Ncell);
Loc = file.param.Loc;
T=ceil(max(s1(1,:)));
Tburn =1000;



%%% check if we have enough neurons to do resampling
FR = hist(s1(2,s1(1,s1(2,:)>0)>Tburn),1:N)/(T-Tburn)*1e3;
idx = find(FR>1);

eligible = zeros(1,Npop);
Nc = zeros(1,Npop);
for jj= 1:Npop
    eligible(1,jj) = sum(idx > Nsum(jj) & idx<= Nsum(jj+1) ); %% count num of eligible neurons to sample 
    if eligible(1,jj) >= ceil(samplesize*(resampleCounts/2))
        %%% less than total eligible, but more than half needed = samplesize
        Nc(jj) = samplesize;
    elseif eligible(1,jj) >= samplesize
        Nc(jj) = eligible(1,jj); %%% less than half needed but more than samepsize
    else eligible(1,jj) < samplesize
        Nc(jj) = 0; %%less than sampsize = 0 eligible to sample
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculations for coherence based on distance
%%% reshuffling the sampled populations
structure = cell(1,resampleCounts);
for rr = 1:resampleCounts
    structure{1,rr} = sprintf('CohCalc_%s',num2str(rr));
end

s1_tmp=s1(:,s1(1,:)>Tburn);
s1_tmp(1,:)=s1_tmp(1,:)-Tburn;

%%%%


for resamp = 1:resampleCounts
    resamp
    tmpstruc = structure{1,resamp};
    Ic=[];
    for pop = 1:Npop  %% grab randam sample of Nc amount
        if nnz(idx<=Nsum(pop+1) &idx>Nsum(pop))>Nc(pop)
            Ic = [Ic, randsample(idx(idx<=Nsum(pop+1) &idx>Nsum(pop)),Nc(pop))];
        else
            Nc(pop) = nnz(idx<=Nsum(pop+1) &idx>Nsum(pop));
            Ic = [Ic,idx(idx<=Nsum(pop+1) &idx>Nsum(pop))];
        end
    end


    Nc_sum = [0 cumsum(Nc)];

    Dx = pdist2(Loc(Ic,1),Loc(Ic,1),'euclidean');
    Dy = pdist2(Loc(Ic,2),Loc(Ic,2),'euclidean');
    Dx(Dx>0.5)=1-Dx(Dx>0.5);
    Dy(Dy>0.5)=1-Dy(Dy>0.5);
    D = sqrt(Dx.^2 + Dy.^2);

    %Ne=param.Ncell(1);

    T=ceil(max(s1(1,:)));
   %%% Tburn=2000;
    dt=1;
    time=1:dt:T;
    % Tw=200;p=100;
    Tw=1000;p=500; %%% window is 1 sec, 1/2 overlap with p

    re=spktime2count(s1_tmp,Ic,dt,(T-Tburn)/dt,0)*1e3/dt; %%% count spk times
    re = re - repmat(mean(re,2),[1,size(re,2)]); %%% subtract mean to reduce noise

    %tic
    L = size(re,2);
    w=ones(Tw,1);
    X=zeros(Tw,ceil(L/(Tw-p)),sum(Nc));
    for mm=1:sum(Nc)
        temp=buffer(re(mm,:),Tw,p);
        X(:,:,mm) = temp.*repmat(w,[1,size(temp,2)]);
    end
    if ceil(L/(Tw-p))*(Tw-p)>L
        X = X(:,2:end-1,:);
    else
        X = X(:,2:end,:);
    end
    Nseg = size(X,2);
    % X=X-repmat(mean(mean(X,1),2),[Tw,Nseg,1]);
    Fx=fft(X);

    %clear X re

    Sx=(dt^2)*1e-3/sum(w.^2)*Fx.*conj(Fx);
    Sx=Sx(1:Tw/2+1,:,:);
    Sx=permute(sum(abs(Sx),2)/(Nseg-1),[1 3 2]); % power spectrum density

    dmax=0.5;
    dd=0.025;
    %%% data storage
    Cd_data = cell(Npop,Npop);
    Cohrxx_mean_data = cell(Npop,Npop);
    Cohrxx_d_data = cell(Npop,Npop);
    Npair_data = cell(Npop,Npop);
    %%%
    D=triu(D,1);
    for jj = 1:Npop %%% row
        for ii = jj:Npop %%% col
            [Cd, Cohrxx_mean, Cohrxx_d, Npair, Nd daxis]= coh_dist(D((1+Nc_sum(jj)):Nc_sum(jj+1), (1+Nc_sum(ii)):Nc_sum(ii+1)), ...
                Sx(:,(1+Nc_sum(jj)):Nc_sum(jj+1)),Sx(:,(1+Nc_sum(ii)):Nc_sum(ii+1)), ...
                Fx(:,:,(1+Nc_sum(jj)):Nc_sum(jj+1)),Fx(:,:,(1+Nc_sum(ii)):Nc_sum(ii+1)), ...
                Tw, Nseg, w, dd,dt, dmax);

            Cd_data{jj,ii} = real(Cd); %% 4x4 cell of 501x20 size
            Cohrxx_mean_data{jj,ii} = Cohrxx_mean; %% 4x4 cell of 501x1 size
            Cohrxx_d_data{jj,ii} = Cohrxx_d;  %% 4x4 cell of 501x20 size
            Npair_data{jj,ii} = Npair; %% %% 4x4 cell of 20x1 size

            %clear Cd Cohrxx_mean Cohrxx_d Npair;
        end
    end


    Fs=1e3/dt;
    f_range=Fs*(0:(Tw/2))/Tw;

   

    CohPwrSpec_Calc.(tmpstruc).Cd_data = Cd_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x daxis 500x20
    CohPwrSpec_Calc.(tmpstruc).Cohrxx_mean_data = Cohrxx_mean_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x 500 x1
    CohPwrSpec_Calc.(tmpstruc).Cohrxx_d_data = Cohrxx_d_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x daxis 


   % clear Cd_data Cohrxx_mean_data Cohrxx_d_data D


end
CohPwrSpec_Calc.eligibleToSamp = eligible;
CohPwrSpec_Calc.Npair_data = Npair_data; %% upper tri cell Npop x Npop, 1 x daxis
CohPwrSpec_Calc.Freq = f_range; %% frequency 0 - 500
CohPwrSpec_Calc.fs = Fs; %% sampling raate
CohPwrSpec_Calc.daxis = daxis; % vec: 1x20
CohPwrSpec_Calc.Twoverlap = p;
CohPwrSpec_Calc.Tw = Tw;
CohPwrSpec_Calc.Nc = Nc;

save(matlabdata_filename,'-append','CohPwrSpec_Calc');


%%%%%%%%%
    function [cd, cohrxx_mean, cohrxx_d, npair, Nd, daxis]=coh_dist(Dmat,s_x1, s_x2, ffx1,ffx2, tw, nseg, window, dd,dt, ddmax)

        daxis=0:dd:ddmax; daxis(1)=10*eps;

        [nnd, ind_d] = histc(Dmat,daxis);  %n: Nd x Nneuron
        nnd=sum(nnd(1:end-1,:),2);
        Nd=length(nnd);
        daxis=daxis(1:end-1)+dd/2;
        s_xx_d=zeros(tw/2+1,Nd);
        cohrxx_d=zeros(tw/2+1,Nd);

        cd=zeros(tw/2+1,Nd);
        npair = zeros(Nd,1);

        for hh=1:Nd
            [row, col]=find(ind_d==hh);
            npair(hh) = length(row);
            if isempty(row)==0 %%% find subset of inputs of pairs at specific distances
                s_xx=(dt^2)*1e-3/sum(window.^2)*ffx1(:,:,row).*conj(ffx2(:,:,col));
                s_xx=permute(sum(s_xx(1:tw/2+1,:,:),2)/(nseg-1),[1 3 2]);
                cohrxx=abs(s_xx).^2 ./ (s_x1(:,row) .* s_x2(:,col)); %% coh in two signals

                cohrxx(isnan(cohrxx))=0;
                s_xx_d(:,hh)=mean(s_xx,2); % pwrspc at prwise dist = hh
                cohrxx_d(:,hh)=mean(cohrxx,2); % 501 x 20;coherence of all freq at prwise dist = hh

                tmp = (s_xx) ./ sqrt((s_x1(:,row) .* s_x2(:,col))); %% 
                tmp(isnan(tmp))=0;
                cd(:,hh)=mean(tmp,2); %%% 501 x 20 % *num pairs/sum(npairs) then needs to be averaged across distances and 
            end
        end
        %toc
        % Sx=mean(Sx,2);
        cohrxx_mean = cohrxx_d*npair/sum(npair); % mean coherence across all distances %% 501 x 1
    end
%%%%%%%%


end





%%%
% re2=re';
%
% tic
% count =0 ;
% freq=0:500;
% Cxy =zeros(size(freq));
% for ii = 1:(Nc(1)-1)
%        tmp = mscohere(re2(:,ii),re2(:,ii+1:Nc(1)),ones(Tw/dt,1),p/dt,freq,1e3/dt);
%        Cxy = Cxy + sum(tmp,1);
%        count = count+size(tmp,1);
% end
% Cxy = Cxy/count;
% toc
%
% figure
% plot(freq,Cxy)
% xlim([0 100])
