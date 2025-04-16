function [syndata] = synchronydata(T,Tburn,neurons,s1,Nsum,Ncell)
%rng(1); %random seed setting to check neurons with same ids
%for testing
%Npop =4;
% e_pop = sort(randi([Nsum(1)+1 Nsum(2)],1,neurons));
% pv_pop = sort(randi([Nsum(2)+1 Nsum(3)],1,neurons));
% som_pop = sort(randi([Nsum(3)+1 Nsum(4)],1,neurons));
% vip_pop = sort(randi([Nsum(4)+1 Nsum(end)],1,neurons));
% %all = [e_pop(1:neurons/2) pv_pop(1:neurons/2) som_pop(1:neurons/2) vip_pop(1:neurons/2)];
sampleids = zeros(4, neurons);

for i = 1:length(Ncell)
    sampleids(i,:) = sort(randi([Nsum(i)+1 Nsum(i+1)],1,neurons));
end

syndata = synchronymeasure(s1,sampleids,T,Tburn, Ncell, Nsum);


%chi_E = synchronymeasure(s1,e_pop,T,Tburn);
%chi_PV = synchronymeasure(s1,pv_pop,T,Tburn);
%chi_SOM = synchronymeasure(s1,som_pop,T,Tburn);
%chi_VIP = synchronymeasure(s1,vip_pop,T,Tburn);
%chi_all = synchronymeasure(s1,all,T,Tburn);

%syndata = [chi_E,chi_PV,chi_SOM,chi_VIP];
