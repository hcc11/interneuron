function [fano_data, fano] = fanofactor(Var,SC,Nc)

Nc_sum = [0 cumsum(Nc)];
%diff = Nc_sum(1,2:end) - Nc;
fano = Var./SC';

fano_data = zeros(1,length(Nc));
for i = 1:length(Nc)
    if Nc(i) ~=0
        fano_data(i) = round(mean(fano(Nc_sum(i)+1:Nc_sum(i+1),1)),3);
    end
end

fano_data = fixnan(fano_data);
% end
% 
% epop_fanomean = round(mean(fano),3);
% pvpop_fanomean = round(mean(pvpop_fano),3);
% sompop_fanomean = round(mean(sompop_fano),3);
% vippop_fanomean = round(mean(vippop_fano),3);
% fano_data =[epop_fanomean, pvpop_fanomean, sompop_fanomean, vippop_fanomean];
