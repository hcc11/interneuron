function [Wrr_new, Wrx_new] = Wremap(Wrr, Wrx, Kout, Kx, Ncell, Nx)

rec_outs = sum(Kout);

Nsum = [0 cumsum(Ncell)];


%%% add second column for the neuron ID providing the current (pre syn,
%%% column 1 is post syn)
rows = 1;
for j =1:length(Ncell)
    for i = Nsum(j)+1: Nsum(j+1)
        count = 1;
        while count <= rec_outs(j)
            Wrr(rows,2)=i;
            count = count +1;
            rows = rows +1;
        end
    end
end
Wrr_new = Wrr;

rows = 1;
for j =1:Nx
    count = 1;
    while count <= sum(Kx)
        Wrx(rows,2)=j;
        count = count +1;
        rows = rows +1;
    end
end
Wrx_new = Wrx;

end