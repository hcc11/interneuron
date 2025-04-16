%%%% function  to calculate synchony measure across network

function [chi] = synchronymeasure(s1,sampleids,T,Tburn,Ncell,Nsum)

%%%% call sdf for each neuron and store it
%sdfs = []; % generate empty cell to fill with spk densitoes for each subneuron
%count = 1;
%tic


%%%%%% this is really slow because it looks at one neuron at a time then
%%%%%% takes the average
%for i=subneurons %%% should be vector value for subneurons
%    [sdfs{count},tt] = spkdensity(s1,i,T,Tburn); %store sdfs and time vector
%    count = count + 1;

%end
%toc

%tic
%%% average sdf
%avg_sdf = zeros(1, length(tt)); %generate emptey vector to full

%for ii=1:length(sdfs)
%    avg_sdf = sdfs{ii} + avg_sdf; % add all sdf's together
%end
%avg_sdf = 1/length(subneurons).*avg_sdf; %average over all subneurons


%time=0:1:T; %% in milliseconds
dt=0.05;  %same dt as sim
time=0:dt:T; %% in milliseconds
re = zeros(length(time),length(Ncell));
re_individ = cell(length(Ncell),1);


for p=1:length(Ncell)
    re_individ{p} = zeros(length(time),size(sampleids,2));
    for i=1:size(sampleids,2)
        re_individ{p}(:,i)=hist(s1(1,s1(2,:)==sampleids(p,i)),time)*1e3; %Hz, depends on dt
    end
    re (:,p) = mean(re_individ{p},2); %average rate across the samples
end


Tw = 200;
re_smoothed=imfilter(re(Tburn/dt:end,:),ones(Tw,1)/Tw); % average rates (average over sample)
re_smoothed_individ = cell(length(Ncell),1); %individual rates
for p=1:length(Ncell)
    re_smoothed_individ{p} = imfilter(re_individ{p}(Tburn/dt:end,:),ones(Tw,1)/Tw);
end

var_individ = zeros(length(Ncell),size(sampleids,2));
chi = zeros(1,length(Ncell));
var_acrosspop = var(re_smoothed,1,1);
for i = 1:length(Ncell) 
    var_individ(i,:) = var(re_smoothed_individ{i},1,1);
    meanindivid = mean(var_individ(i,:),2);
    if meanindivid == 0
        chi(i) = 0;
    else
        chi(i) = sqrt(var_acrosspop(1,i)/meanindivid);
    end
end





%%%% sigmavar of individual neurons
%svar = {}; % create emtpy cell to store each varience for sdf
%avgii_svar = 0; 
%for j = 1:length(sdfs) %for each sfts
%    svar{j} = sigmavar(sdfs{j}, tt);  %calculate the sigma variance 
%    avgii_svar = avgii_svar + svar{j}; % add all sigma varience for individual sdf together
%end

%avgii_svar = 1/(length(sampleids)).*avgii_svar; %average all sigma variance for each sdf

%%%% average sigmavar: sigmavar of all neurons in pop
%total_svar = sigmavar(avg_sdf,tt); %varience over entire set of subneurons

%%% synchony measurement
%chi = sqrt(total_svar/avgii_svar); 
%toc


end