function [IdsX_final_cell] = Wframe_PresynId(Wrx, Kxout, Nx, sampids)
IdsX_final_cell=cell(1,1);
IdsX_final = sortrows(Wrx(ismember(Wrx(:,1),sampids),:),[1,2]); %put input pop connections in individual cell
IdsX_final = double(IdsX_final); % necessary for hist call in presynspks.m
IdsX_final_cell{1,1} = IdsX_final; 
clear IdsX_final