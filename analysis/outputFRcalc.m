function [FR] = outputFRcalc(sampids_pop, T, Tburn, s1)

s1burned = s1(2,s1(1,:) > Tburn); % limit spk times preemptively 
s1burned = sort(s1burned(ismember(s1burned,sampids_pop)));  % find and count spks
FR = arrayfun(@(x) 1000*nnz(ismember(s1burned,x)),sampids_pop) ./ (T-Tburn); %calc avg firing rate

% varFR = zeros(1,numel(sampids_pop));
% s1sm = s1(:,s1(1,:) > Tburn);
% Tw = 50;
% for ss = 1:numel(sampids_pop)
%     tmp=hist(s1sm(1,s1sm(2,:)==sampids_pop(ss)),Tburn:T)*1e3;
%     tmp2 = imfilter(tmp,ones(Tw,1)/Tw);
%     varFR(1,ss) = var(tmp2);
%     clear tmp tmp2
% end

% step = 1;
% time=0:step:T;
%    % time=t1:step:t2;plot(tmp2)
% re = zeros(length(time),Npop); %data allocated space
% for p=1:Npop
%     re(:,p)=hist(srr(1,srr(2,:)<=Nsum(p+1)&srr(2,:)>Nsum(p)),time)/Ncell(p)*1e3/step; % firing rate for each pop
% end
% Tw = 50;
% re_smoothed=imfilter(re,ones(Tw,1)/Tw); %filter firing rate
% nuSimPoisson = mean(re,1);

end

% Original from 2 March 2023
% function [FR] = outputFRcalc(sampids_pop, T, Tburn, s1)
% 
% s1burned = s1(2,s1(1,:) > Tburn); % limit spk times preemptively 
% s1burned = sort(s1burned(ismember(s1burned,sampids_pop)));  % find and count spks
% FR = arrayfun(@(x) nnz(ismember(s1burned,x)),sampids_pop) ./ (T-Tburn); %calc avg firing rate
% 
% end
