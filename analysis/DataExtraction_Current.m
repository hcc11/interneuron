%%%%%% pull data from trial runs and calculations into two matlab files:
%%%%%%%
%%%%%%%         uses 'paramVal1' in each file as param values. 
function DataExtraction_Current(data_folder, filename, like_filename)
%function DataExtraction(option, data_folder, filename, like_filename, param_rowm, param_coln)
%%%%%%% saves a file as ['AllCurrent_Avgs_', filename] in data_folder 


%if ~isfield(option, 'dataextract') option.dataextract=0; end
%if option.dataextract == 1
% % Initial setup: load files
%like_newfilename = erase(like_filename, 'ID1')
like_newfilename = like_filename;

original_files = dir([data_folder like_newfilename])
file_num = length(original_files)
Npop = 4;


% matlabdata_filename_individ = sprintf("%s%s%s",data_folder,'AllData_Individ_', filename);
matlabdata_filename_avgs = sprintf("%s%s%s",data_folder,'AllCurrent_Avgs_', filename);
%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%
%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%

tic;
%if option.acrosstrialcalcs == 1
disp('Current Data Grab')
count = 1
for ff = 1:file_num

    ff
    % % load .mat simulation file
    file = load([data_folder original_files(ff).name]);
     original_files(ff).name

    if contains(original_files(ff).name, 'ID1') == 1
        IndividData.currentCalc{count} = file.currentCalc;
       %%%% static input
        IndividData.paramval(count) = file.paramVal1; 
        IndividData.nuSim{count} = file.nuSim;
        count = count + 1;
        if ff == 1
            IndividData.param = file.param;
        end
    end

    %%%%%%%% current data %%%%%%
    currentdata.totalMean{ff} = file.currentCalc.total_mean;
    currentdata.totalVar{ff} = file.currentCalc.total_var;

    currentdata.totalEmean{ff} = file.currentCalc.E_mean;
    currentdata.totalEvar{ff} = file.currentCalc.E_var;

    currentdata.totalImean{ff} = file.currentCalc.I_mean;
    currentdata.totalIvar{ff} = file.currentCalc.I_var;

    currentdata.nuSim{ff} = file.nuSim; %average firing rates recorded
   
    %%%% NOTE: This needs to be changed to the correct parameter that
    %%%% changes across trials
    data.alpha(ff) = file.paramVal1;
    %%%% static == {param_rowm, param_coln} = {12,2}
    %%%%heterogen === {param_rowm, param_coln} = {13,2}(popwithinput)
   %  data.alpha(ff) = file.ParamChange{param_rowm, param_coln};
   % data.alpha(ff) = file.ParamChange{13, param_coln}(popwithinput);

end
%save(matlabdata_filename_individ, 'CorrCurrentData_eachtrial','CorrFRData_eachttrial','-v7.3');
save(matlabdata_filename_avgs, 'IndividData','-v7.3');


%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%
%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%


%%%% save multicalcs
Npop = 4;
samplesize = 500;
unique_param = unique(data.alpha); %finds unique values of changing parameter

nuSim_mean_Current = cell(1,length(unique_param));
output_data = cell(4,length(unique_param)); %row 1: FR mean output;

TotalMean_Current = cell(Npop,length(unique_param));
TotalVar_Current = cell(Npop,length(unique_param));

TotalMean_ECurrent = cell(Npop,length(unique_param));
TotalVar_ECurrent = cell(Npop,length(unique_param));
TotalMean_ICurrent = cell(Npop,length(unique_param));
TotalVar_ICurrent = cell(Npop,length(unique_param));

for ss = 1:length(unique_param)
    matches = data.alpha == unique_param(ss); % Grab data.alpha locations that match specific alpha
    match = 0; % Index to see the matches (running count)
    for kk = 1:length(matches) % Run through all possible matches to look for non zero matches
        if matches(kk) == 1 % If non zero match, add to file_match struct
            match = match + 1; % Increase index
            
            file_match.nuSim_mean{match} = currentdata.nuSim{kk};

            file_match.totalmean{match} = currentdata.totalMean{kk};
            file_match.totalvar{match} = currentdata.totalVar{kk};

            file_match.Emean{match} = currentdata.totalEmean{kk};
            file_match.Evar{match} = currentdata.totalEvar{kk};
            file_match.Imean{match} = currentdata.totalImean{kk};
            file_match.Ivar{match} = currentdata.totalIvar{kk};

        end
    end
    % At the point in the code: all the files that correspond to a single changing
    % parameter value (i.e different trials across the same
    % condition) are collected.  Next, calculate the averages
    % across all trials for one condition

    Var_tmp = zeros(length(file_match.totalmean),4);
    Mean_tmp = zeros(length(file_match.totalmean),4);

    Evar_tmp = zeros(length(file_match.totalmean),4);
    Emean_tmp = zeros(length(file_match.totalmean),4);

    Ivar_tmp = zeros(length(file_match.totalmean),4);
    Imean_tmp = zeros(length(file_match.totalmean),4);

    nuSim_mean_tmp = zeros(length(file_match.totalmean),4);

    for mm = 1:length(file_match.totalmean)
        for nn = 1:4
            Var_tmp(mm,nn) = mean(file_match.totalvar{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));
            Mean_tmp(mm,nn) = mean(file_match.totalmean{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));

            Evar_tmp(mm,nn) = mean(file_match.Evar{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));
            Emean_tmp(mm,nn) = mean(file_match.Emean{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));

            Ivar_tmp(mm,nn) = mean(file_match.Ivar{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));
            Imean_tmp(mm,nn) = mean(file_match.Imean{1,mm}((nn-1)*samplesize+1:nn*samplesize,1));
        end
        nuSim_mean_tmp(mm,:) = file_match.nuSim_mean{mm};
    end
    nuSim_mean(ss,:) = mean(nuSim_mean_tmp,1);
    Mean_Current(ss,:) = mean(Mean_tmp,1);
    Var_Current(ss,:) = mean(Var_tmp,1);
    Emean_Current(ss,:) = mean(Emean_tmp,1);
    Evar_Current(ss,:) = mean(Evar_tmp,1);
    Imean_Current(ss,:) = mean(Imean_tmp,1);
    Ivar_Current(ss,:) = mean(Ivar_tmp,1);
    paramVal(ss,:) = unique_param(ss);

end

corrdata_current.paramVal = paramVal;
corrdata_current.nuSim_avg = nuSim_mean;
corrdata_current.TotalMean_Current = Mean_Current;
corrdata_current.TotalVar_Current = Var_Current;
corrdata_current.Emean_Current = Emean_Current;
corrdata_current.Evar_Current = Evar_Current;
corrdata_current.Imean_Current = Imean_Current;
corrdata_current.Ivar_Current = Ivar_Current;

save(matlabdata_filename_avgs, 'corrdata_current', '-append');

end

