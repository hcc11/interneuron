function [totalinputcurrent EI_totalinput output_Each output_All] = currentandspks_calc(presynids, presynXids, samppopids, ...
    samppoptype, s1, sx, scaledweights, tauds, taurs, time, dt, samplesize, Iapp, Nsum, Nx, ouMean, ouSigma)

% tau_r = taurs(1);
step=20;

filterAll= cell(1,5);

for ii = 1:5
    tt = -10*tauds(ii):dt:10*tauds(ii); %% filter adapts based on tau_decay
    filter =(exp(-tt./tauds(ii))-exp(-tt./taurs(ii))).*(tt>=0);
    filterAll{1,ii} = filter/sum(filter)/dt;
    clear tt filter
end


%%%% matrix allocations
totalinputcurrent = zeros(samplesize, length(time(1:step:end)));
EI_totalinput = zeros(2, length(time(1:step:end))); % excit(I>0) current row 1, inhib (I<0) current row 2
output_Each = zeros(samplesize,2);

for jj = 1:samplesize
    nodupe = unique(presynids{jj,1}(:,2));
    if ~isempty(nodupe) == 1

        %%% now we count and calc duplicate connections
        dupecount = [nodupe,(-1 + hist(presynids{jj,1}(:,2),nodupe))'];
        singledupe = dupecount(dupecount(:,2)==1,1);
        multidupe = dupecount(dupecount(:,2)>1,:);

        for prepop = 1:length(Nsum)-1
            if isempty(nodupe(nodupe(:,1) <= Nsum(prepop+1) & nodupe(:,1) >= Nsum(prepop)+1))==0
                spks_tmp1 = s1(1,ismember(s1(2,:), nodupe(nodupe(:,1) <= Nsum(prepop+1) & nodupe(:,1) >= Nsum(prepop)+1)));
                spks_single = s1(1,ismember(s1(2,:),singledupe(singledupe(:,1) <= Nsum(prepop+1) & ...
                    singledupe(:,1) >= Nsum(prepop)+1)));
                multispks = [];
                multidupe_prepop = multidupe(multidupe(:,1)<= Nsum(prepop+1) &multidupe(:,1)>= Nsum(prepop)+1,:);
                if ~isempty(multidupe_prepop) == 1

                    for kk = 1:size(multidupe_prepop,1) % multi connection dupes
                        multispks = [multispks repmat(s1(1,s1(2,:)== multidupe_prepop(kk,1)),1,multidupe_prepop(kk,2))]; % find spk times
                    end
                    %                 spks_tmp = [spks_tmp multispks];
                end
                spks_tmp = [spks_tmp1 spks_single multispks];
                %             tmp_current = zeros(1,length(time(1:step:end)));
                if isempty(spks_tmp) == 0
                    tmp_current = hist(spks_tmp(spks_tmp<time(end) & spks_tmp>time(1)),time); % ms
                    tmp_current = imfilter(tmp_current, filterAll{1,prepop}, 'conv');
                    tmp_current = tmp_current(1:step:end).*scaledweights(prepop);
                    totalinputcurrent(jj,:) =  totalinputcurrent(jj,:) + tmp_current;
                    if prepop ==1
                        EI_totalinput(1,:) =  EI_totalinput(1,:) + tmp_current;

                    elseif prepop > 1
                        EI_totalinput(2,:) =  EI_totalinput(2,:) + tmp_current;
                    end
                    %  clear tmp_current
                end
            end
        end
    end
    %step = 20; % save every 1 sec of input
    %%% tmp_current = tmp_current save subset of current with step =
    %%% 1
    %%% change to matrix becuase we dont need to save prepop just save total
    %%%  inputcurrent{1,prePop}(jj,:)= tmp_current.*scaledweights(prePop);


    outspks_tmp = s1(1,ismember(s1(2,:), samppopids(jj)));
    tmp = hist(outspks_tmp(outspks_tmp<time(end) & outspks_tmp>time(1)),time); % ms
    output_Each(jj,1) = 1000*nnz(tmp)/time(end);
    output_Each(jj,2) = scaledweights(samppoptype)*output_Each(jj,1);

end
clear nodupe singledupe multidupe spks_tmp tmp_current outputspks_tmp tmp


outall_tmp = s1(1,ismember(s1(2,:), samppopids));
tmp = hist(outall_tmp(outall_tmp<time(end) & outall_tmp>time(1)),time); % ms
output_All(1,jj) = 1000*nnz(tmp)/time(end);
output_All(2,jj) = scaledweights(samppoptype)*output_All(1,jj);




%  totalinputcurrent(jj,:) = tmp_current;

if ~isempty(sx) == 1
    prepop = 5;

    for jj = 1:samplesize

        nodupe = unique(presynXids{jj,1}(:,2));
        if isempty(nodupe)==0
            %%% now we count and calc duplicate connections
            dupecount = [nodupe,(-1 + hist(presynXids{jj,1}(:,2),nodupe))'];
            singledupe = dupecount(dupecount(:,2)==1,1);
            multidupe = dupecount(dupecount(:,2)>1,:);

            spks_tmp1 = sx(1,ismember(sx(2,:), nodupe(:,1)));
            spks_single = sx(1,ismember(sx(2,:),singledupe(:,1)));

            multispks = [];
            if ~isempty(multidupe) == 1
                for kk = 1:size(multidupe,1) % multi connection dupes
                    multispks = [multispks repmat(sx(1,sx(2,:)== multidupe(kk,1)),1,multidupe(kk,2))];
                end

            end
            spks_tmp = [spks_tmp1 spks_single multispks];

            if isempty(spks_tmp) == 0

                %    tt = -10*tauds(prepop):dt:10*tauds(prepop); %% filter adapts based on tau_decay
                %    filter =(exp(-tt./tauds(prepop))-exp(-tt./tau_r)).*(tt>=0);
                %    filter = filter/sum(filter)/dt;

                tmp_current = hist(spks_tmp(spks_tmp<time(end) & spks_tmp>time(1)),time); % ms
                tmp_current = imfilter(tmp_current, filterAll{1,prepop}, 'conv');
                tmp_current = tmp_current(1:step:end).*scaledweights(prepop);
                totalinputcurrent(jj,:) = totalinputcurrent(jj,:) + tmp_current;
                EI_totalinput(1,:) =  EI_totalinput(1,:) + tmp_current; %% x network is excitatory
            end

        end

        %step = 20; % save every 1 sec of input
        %%% tmp_current = tmp_current save subset of current with step =
        %%% 1
        %%% change to matrix becuase we dont need to save prepop just save total
        %%%  inputcurrent{1,prePop}(jj,:)= tmp_current.*scaledweights(prePop);

    end


end

totalinputcurrent = totalinputcurrent + Iapp;
if Iapp > 0
    EI_totalinput(1,:) =  EI_totalinput(1,:) + Iapp;
elseif Iapp <= 0
    EI_totalinput(2,:) =  EI_totalinput(2,:) + Iapp;
end

if ~isempty(ouMean)==1 & ~isempty(ouSigma)==1

    for jj = 1:samplesize

        %ouMean, ouSigma
        tau = 5;    % Time constant in ms
        
        xx = zeros(1, length(time(1:step:end)));
        xx(1) = randn; %random ic
        % Generate the OU process
        for kk = 2:length(time(1:step:end))
            dW = sqrt(dt)*randn;
            xx(kk) = xx(kk-1)+(-1/tau)*(xx(kk-1) - ouMean)*dt + (ouSigma/tau)*dW;
        end

      %  mean(xx), var(xx)

      %  mean(totalinputcurrent(jj,:))
      %  var(totalinputcurrent(jj,:))
        
        totalinputcurrent(jj,:) = totalinputcurrent(jj,:) + xx;

      %  mean(totalinputcurrent(jj,:))
      %  var(totalinputcurrent(jj,:))
        clear xx


    % if ouMean > 0 
    %     EI_totalinput(1,:) =  EI_totalinput(1,:);
    % elseif ouMean < 0
    %     EI_totalinput(2,:) =  EI_totalinput(2,:);
    % else ouMean = 0
    %     %%% I need to add what type of excitatory or inhibitory current
    %     %%% this would be or EI splitting of current - but not critical
        %%% right now
    %end

    end


end




EI_totalinput = EI_totalinput/samplesize;

%%% estimated total average current from sampled pop (ie output
%%% current)

end

