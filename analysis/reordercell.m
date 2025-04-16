function reordered = reordercell(presyn,sampids,samplesize)

%if size(presyn,2) == 4
    reordered = cell(samplesize,1);
    for nn = 1:samplesize
        tmp = [];
        for gg = 1:size(presyn,2)
             tmp = [tmp; presyn{1,gg}((presyn{1,gg}(:,1) == sampids(nn)),:)];
        end
      reordered{nn,1} = tmp;
    end
%end
