function [maxfreq] = peakfreq(psdx_all)
% function [e_freqratio, pv_freqratio, som_freqratio, vip_freqratio] = peakfreq_ratio(psdx_all)
% output: ratio of max peak freuquency to 0th frequency. ratio close to 0 =
% nonosicllatory, if ratio >=0.2 oscillatory population dynamics (for sure???)

%%%%% peak power ratio with 0 frequency ratio
%%% E pop
e_0freq = psdx_all(1,1); %grab 0th freq peak
[epks,elocs] = findpeaks(psdx_all(:,1)); %find all peaks and x coord in powerspectrum
%e_pkcoord = [epks elocs]; %combine to new vec (not needed)
e_maxfreq = max(epks); %find max peak value

%e_freqratio = e_maxfreq/e_0freq; %ratio of max peak over 0 freq


%%% PV pop
pv_0freq = psdx_all(1,2);
[pvpks,pvlocs] = findpeaks(psdx_all(:,2));
pv_pkcoord = [pvpks pvlocs];
pv_maxfreq = max(pvpks);

%pv_freqratio = pv_maxfreq/pv_0freq;


%%% SOM pop
som_0freq = psdx_all(1,3);
[sompks,somlocs] = findpeaks(psdx_all(:,3));
som_pkcoord = [sompks somlocs];
som_maxfreq = max(sompks);

%som_freqratio = som_maxfreq/som_0freq;


%%% VIP pop
vip_0freq = psdx_all(1,4);
[vippks,viplocs] = findpeaks(psdx_all(:,4));
vip_pkcoord = [vippks viplocs];
vip_maxfreq = max(vippks);

%vip_freqratio = vip_maxfreq/vip_0freq;

if isempty(e_maxfreq)
    e_maxfreq = 0;
end
if isempty(pv_maxfreq)
    pv_maxfreq = 0;
end
if isempty(som_maxfreq)
    som_maxfreq = 0;
end
if isempty(vip_maxfreq)
    vip_maxfreq = 0;
end

maxfreq = [e_maxfreq, pv_maxfreq, som_maxfreq, vip_maxfreq];

end