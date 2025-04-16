 function [presynIDs sampids] = Wframe_PresynId(Wrr, Kout, Ncell, poptype, samplesize, sampids)

%Nsum = cumsum([0 Ncell]);
%sampids = 31461; %% testing recorded ID
%sampids = sort([31461, 37833, 1782, 14399]); testing recorded set
%sampids = sort(Nsum(poptype)+ randperm(Ncell(poptype),samplesize)); %unique sample of IDs
n=sum(Kout).*Ncell; %% total number of connections
bounds = [0 cumsum(n)]; %% boudaries for Wrr matrix corresponding to the first index
                      %%% of where the input type changes to a different
                      %%% pop type

presynIDs = cell(1,length(Ncell));
for j = 1:length(Ncell) % for each input pop
    membership = ismember(Wrr(bounds(j)+1:bounds(j+1),1),sampids);
    subWrr = Wrr(bounds(j)+1:bounds(j+1),:); % sub matrix to reindex all Ids (easier to grab from ismemeber)
    presynIDs{1,j} = double(sortrows(subWrr(membership,:),[1,2])); %put input pop connections in individual cell
                % % double needed for hist call later in presynspks.m
end

