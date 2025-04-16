function Calc_corrFR(data_folder, like_filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('MultiPopSpatialNet Calc')

original_files = dir(like_filename);
file = load([data_folder original_files.name]);
matlabdata_filename = like_filename;
allnames = split(file.ParamChange{1,2},'/');
datafilename = char(allnames(end));

s1 = file.s1;

Npop = file.ParamChange{8,2};
N = file.param.N;
Nx = file.param.Nx;
Nc = file.ParamChange{4,2};
Ncell = file.param.Ncell;
Nsum = [0, cumsum(Ncell)];
Loc = file.param.Loc;
nuSim = file.nuSim;

T = file.param.T;
Tburn =2000;
Jr = file.param.Jr;
Jx = file.param.Jx;

param = file.paramVal1;

binsize=0.5;
samplesize = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlation calculation for single simulation
% if option.SingleCalc == 1

[C_d,COV_d,Cbar,COVbar,daxis,Var,SC,Tw,Nc]=corr_d(s1,Ncell,Loc,Nc);
Nc_sum = [0 cumsum(Nc)];

[fanomean, fano] = fanofactor(Var,SC,Nc);

[chi] = synchronydata(T,Tburn,samplesize,s1,Nsum, Ncell);

[psdx, freq, re_unsmoothed] = power_spectrum(s1,T,Ncell,Npop);

[peaks] = peakfreq(psdx);

corrCalc.C_d = C_d;
corrCalc.COV_d = COV_d;
corrCalc.Cbar = Cbar;
corrCalc.COVbar = COVbar;
corrCalc.daxis = daxis;
corrCalc.Var = Var;
corrCalc.SC = SC;
corrCalc.Tw = Tw;
corrCalc.Nc = Nc;
corrCalc.Nc_sum = Nc_sum;
corrCalc.fanomean = fanomean;
corrCalc.fano = fano;
corrCalc.psdx = psdx;
corrCalc.freq = freq;
corrCalc.chi = chi;
corrCalc.peakdata = peaks;
corrCalc.reUnsmoothed = re_unsmoothed;
corrCalc.param = param;
corrCalc.Loc = file.param.Loc;


save(matlabdata_filename,'corrCalc','-append');

clear all;
end