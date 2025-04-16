%%%%% calculated frequency spectrum of neurons
function [psdx_all, freq, re]=power_spectrum(s1,T,Ncell,Npop)
%%%% psdx : powerspectrum data
data_structure={};
Tburn=200;
dt = 1;
time = 0:dt:T;

re = zeros(length(time),Npop);
Nsum = [0 cumsum(Ncell)];


%%%%% rate is spks per timebin (based on time = 0:dt:T) therefore spks/ms *1e3 = spks/Hz
%%%% divide but number of cells in each pop since avgerage over all cells
re = zeros(length(time),Npop);
for p=1:Npop
    re(:,p)=hist(s1(1,s1(2,:)<=Nsum(p+1)&s1(2,:)>Nsum(p)),time)*1e3/Ncell(p); %% Hz
end

re = re(Tburn/dt+1:end,:);
%re_smoothed=re(Tw/2-1:end-Tw/2,:);

data_structure{1}= re;  %collect population firing rate data 
xx = zeros(length(data_structure{1}),4);
psdx_all = zeros(ceil(length(xx)/2),4);

fs =1/(dt/1000);% dt/1000 convert ms to sections % dt == 1ms, fs = 1/dt, sampling at 1000Hz
N = length(xx); %number of samples
freq = 0:fs/N:fs/2;
xx = data_structure{1}; %E, PV, SOM, VIP

for jj=1:Npop

    xdft = fft(xx(:,jj));
    xdft = xdft(1:ceil(length(xx)/2));
    psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    psdx_all(:,jj) = psdx;

%     %%%%%% test fft against periodgram
%
%     %%%[p,f] = pspectrum(xx(:,jj),fs); % dont user pspectrum not the same, power (Watts)
%     [p2,f2] = periodogram(xx(:,jj),rectwin(N),N,fs); %% power spec density function (Watts/Hz), same as fft
%     
%     subplot(1,3,1)
%     plot(freq(2:end),psdx(2:end));
%     title('FFT calculated Power Spec');
%     xlim([0,100]);
%     subplot(1,3,2)
%     plot(f2(2:end),p2(2:end))
%     xlim([0,100]);
%     title('Peridogram calculated Power Spec');
%     subplot(1,3,3)
%     hold on
%     plot(f2(2:end),p2(2:end), '--')
%     plot(freq(2:end),psdx(2:end));
%     xlim([0,100]);
%     title('Plotting Overlap Power Spec');
%     %legend({'fft','pspectrum'});
%     hold off
%     %psdx_all(:,jj) = psdx;  %%% hz
end

