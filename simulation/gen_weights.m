function [Wrr,Wrf, K, Kx]=gen_weights(Ncell,Nx,sigmaRX,sigmaRR,Prr,Prx,Loc,X_Loc)
% takes ~280 sec to run for N = 5e4
% Wrr, Wrf  : save postsynaptic index for each neuron
% sort post syn index

% Loc: (x, y) index
% Ncell:  # of cells per population in the recurrent layer
Ncell = reshape(Ncell,1,[]);
Npop = length(Ncell) ;
Nidx = [0 cumsum(Ncell)];

% K(i,j) is the total number of projections a neuron from pop j makes to ALL neurons in pop i
K = ceil(Prr .* repmat(Ncell(:),1,Npop));
Kx = ceil(Prx.*Ncell(:));
%Kx= Prx.*Ncell(:);  % total number of projections in ffwd layer
%Kx
Wrr=zeros(sum(sum(K,1).*Ncell),1,'int32'); % recurrent connections weight
Wrf=zeros(sum(Kx)*Nx,1,'int32'); % feedforward connections weight

% recurrent connections
count = 0;
for j = 1: Npop  % pre- population
    for id_pre = (1:Ncell(j))+Nidx(j) % index for presynaptic neurons
        Dx = Loc(:,1) - Loc(id_pre,1); % distance in x of current rec neuron to presyn neuron
        Dy = Loc(:,2) - Loc(id_pre,2); % distance in y of current rec neuron to presyn neuron
        for i = 1: Npop  % post- population
            if K(i,j)
                id_post = (1:Ncell(i))+Nidx(i);
                if sigmaRR(i,j)<0.4
                    nwrap = 1;
                else
                    nwrap = 3;
                end
                Prob=WrappedGauss1D(Dx(id_post),0,sigmaRR(i,j),nwrap).*WrappedGauss1D(Dy(id_post),0,sigmaRR(i,j),nwrap);
                % nwrap needs to be large for sigma, but takes longer time
                postcell = randsample(id_post,K(i,j),true,Prob); % sample from K neurons from id_post with probability Prob.
                Wrr(count+(1:K(i,j)))=sort(postcell);
                count = count + K(i,j);
            end
        end
    end
end
count,

% feedforward connections
count=0;
for id_pre = 1:Nx % ffwd populations
    Dx = Loc(:,1) - X_Loc(id_pre,1); % distance in x of ffwd neuron connections to presyn neuron
    Dy = Loc(:,2) - X_Loc(id_pre,2); % distance in y of ffwd neuron connections to presyn neuron
    for i = 1: Npop  % post- population
        if Kx(i)
            id_post = (1:Ncell(i))+Nidx(i);
            if sigmaRX(i)<0.4
                nwrap = 1;
            else
                nwrap = 3;
            end
            Prob=WrappedGauss1D(Dx(id_post),0,sigmaRX(i),nwrap).*WrappedGauss1D(Dy(id_post),0,sigmaRX(i),nwrap);
            postcell = randsample(id_post,Kx(i),true,Prob); % sample from K neurons from id_post with probability Prob.
            Wrf(count+(1:Kx(i)) )=sort(postcell);
            count = count + Kx(i);
        end
    end
end
count,
