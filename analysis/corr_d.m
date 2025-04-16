
function [C_d,COV_d,Cbar,COVbar,daxis,Var,SC,Tw,Nc]=corr_d(s1,Ncell,Loc,Nc)
% s1: spike times of neuron population, 2xN
% N: total # of neurons in all pops {E, PV, SOM, VIP} respectively
% Nc: 1x4 # of sampled neurons, e.g. Nc=[500, 500, 500, 500];
% only neuron in the center square [.25, .75]x[.25, .75] are sampled
% Cbar=[C22, C12, C11], mean correlations between populations
% C_d: [20x3] correlation for each distance (daxis)
% COV_d, COVbar: same as (C_d and Cbar) but for covariance
% SC: mean rate (Hz) of the Nc sampled neurons
% Var: variance of the Nc sampled neurons

T=ceil(max(s1(1,:)));
Tburn =500;
N = sum(Ncell); 
Npop = length(Ncell); 
FR =hist(s1(2,s1(1,s1(2,:)>0)>Tburn),1:N)/(T-Tburn)*1e3; 
idx = find(FR>1);  
% idx=find(FR'>1 & Loc(:,1)<0.75 & Loc(:,1)>0.25 & Loc(:,2)<0.75 & Loc(:,2)>0.25);
% idx = reshape(idx,1,[]); 
Nsum=[0 cumsum(Ncell)]; %cumulative sum: list of elements of cumulative sum; Ncell is num of neurons/population type
Ic=[];
for pop = 1:Npop 
    if nnz(idx<=Nsum(pop+1) &idx>Nsum(pop))>Nc(pop)
        Ic = [Ic, randsample(idx(idx<=Nsum(pop+1) &idx>Nsum(pop)),Nc(pop))];
    else
        Nc(pop) = nnz(idx<=Nsum(pop+1) &idx>Nsum(pop));
        Ic = [Ic,idx(idx<=Nsum(pop+1) &idx>Nsum(pop))];  
    end
end

% distance calculations
Dx = pdist2(Loc(Ic,1),Loc(Ic,1),'euclidean');
Dy = pdist2(Loc(Ic,2),Loc(Ic,2),'euclidean');
Dx(Dx>0.5)=1-Dx(Dx>0.5); 
Dy(Dy>0.5)=1-Dy(Dy>0.5); 
D = sqrt(Dx.^2 + Dy.^2);

% D = pdist2(Loc(Ic,:),Loc(Ic,:),'euclidean');

% compute spike counts using sliding window
%Tw=200; % sliding window size
Tw = 100;
%time=0:1:T;
step=1;   
time = 0:step:T;

re=zeros(length(time),sum(Nc));

for mm=1:sum(Nc)
    re(:,mm)=hist(s1(1,Ic(mm)-1/4<s1(2,:) & s1(2,:)<=Ic(mm)+1/4),time);
end

re_s=imfilter(re(Tburn/step+1:end,:),ones(Tw/step,1)); % re_s is filter for re
%re_s=imfilter(re(Tburn+1:end,:),ones(Tw,1)); % re_s is filter for re
% imfilter(A, H, options)- filters A with H filter, results in size(re_s)==size(H)

re_s=re_s(Tw/step-1:end-Tw/step,:); % 


COV=cov(re_s); % covariance of re_s % produces square matrix 
Var=diag(COV); % varience = diag of of Covariance  % vector
SC=mean(re_s,1); % mean of columns of re_s % vector

R=corrcov(COV);

dmax=0.5; 
dd=0.025;
U=triu(ones(size(R)),1);
Cbar = zeros(Npop,Npop); 
COVbar = zeros(Npop,Npop); 
C_d = cell(Npop,Npop); 
COV_d = cell(Npop,Npop); 
Nc_sum = [0; cumsum(Nc(:))];
for ii =1:Npop
    for jj = ii : Npop
        Utemp=zeros(size(U));
        Utemp((1+Nc_sum(ii)):Nc_sum(ii+1), (1+Nc_sum(jj)):Nc_sum(jj+1))=U((1+Nc_sum(ii)):Nc_sum(ii+1), (1+Nc_sum(jj)):Nc_sum(jj+1));  % index for corr. within recurrent layer (E & I)
        [Crr, COVrr, cbar_rr, COVbar_rr, daxis]=corr_dist(R,COV,D,Utemp, dd, dmax);
        Cbar(ii,jj) = cbar_rr;
        C_d{ii,jj}= Crr;
        COVbar(ii,jj) = COVbar_rr;
        COV_d{ii,jj}= COVrr;
        
    end
end
end




%%%%%%%%%
function [c, cov_m,cbar, COVbar,daxis]=corr_dist(R,COV,D,U, dd, dmax)
% sort corr by distance according to index matrix U
% U is of same size as R & D
R=R(U==1);
D=D(U==1);
COV=COV(U==1);
% dmax=0.5;
% dd=0.025;
daxis=0:dd:dmax;
[n, ind] = histc(D,daxis);
n=n(1:end-1); % discard last bin (d>dmax)
c=zeros(length(n),1);
cov_m=zeros(length(n),1);
for k=1:length(n)
    c(k)=mean(R(ind==k));
    cov_m(k)=mean(COV(ind==k));
end
daxis=daxis(1:end-1)+dd/2;  % center pts
cbar=mean(R(D<=dmax)); % average corr.
COVbar=mean(COV(D<=dmax));
end



