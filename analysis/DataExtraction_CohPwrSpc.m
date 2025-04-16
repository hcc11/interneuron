%%%%%% pull data from trial runs and calculations into two matlab files:
%%%%%%%         File (1): Individual calculations for each trial
%%%%%%%         File (2): Trial Averages, cacluations across trials of the
%%%%%%%         same parameters collected with other averages at different
%%%%%%%         param values. 
%%%%%%%         uses 'paramVal1' in each file as param values. 

%%%%%%% resampcount === number of resampled coh calcs written as
%%%%%%% file.CohPwrSpc_Calc.CohCalc_XX where XX = number 1:resampcount
function DataExtraction_CohPwrSpc(data_folder, filename, like_filename, resampcount)
%%%%%%% saves a file as ['AllCohPwr_Avgs_resamp_', filename] in data_folder 

%%%%% resamp is the number of times that the coherence calculation is
%%%%% uniquely done for each file 

%if ~isfield(option, 'dataextract') option.dataextract=0; end
%if option.dataextract == 1
% % Initial setup: load files
%like_newfilename = erase(like_filename, 'ID1')
like_newfilename = like_filename;

original_files = dir([data_folder like_newfilename]);
file_num = length(original_files);


% matlabdata_filename_individ = sprintf("%s%s%s",data_folder,'AllCohPwr_Individ_', filename);
matlabdata_filename_avgs = sprintf("%s%s%s",data_folder,'AllCohPwr_Avgs_resamp_', filename)
%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%
%%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%

tic;
%if option.acrosstrialcalcs == 1
disp('Coherence Pwr Spec')
count = 1 %%% count keeps track of where to put individual data from each trial== ID1, not the same as resampcount
for ff = 1:file_num
    ff
    % % load .mat simulation file
    file = load([data_folder original_files(ff).name]);
    original_files(ff).name
    Npop = 4;
    daxis = file.CohPwrSpec_Calc.daxis;
    freq = file.CohPwrSpec_Calc.Freq;
    

    

    if resampcount>1
        %%%%%
        ResampAvg.daxis = daxis;
        ResampAvg.Freq = freq;    
        
        ResampAvg.Npair_data = file.CohPwrSpec_Calc.Npair_data;
        ResampAvg.fs = file.CohPwrSpec_Calc.fs;
        ResampAvg.Twoverlap = file.CohPwrSpec_Calc.fs;
        ResampAvg.Tw = file.CohPwrSpec_Calc.Tw;
        ResampAvg.Nc = file.CohPwrSpec_Calc.Nc;
        %%%%%
        
        Cd_mean_Internal = cell(Npop,Npop);
        Cohmean_Internal = cell(Npop,Npop);
        Cohd_Internal = cell(Npop,Npop);

        for sampcnt = 1:resampcount
            tmpname = sprintf('CohCalc_%s',num2str(sampcnt));

            tmpdata = file.CohPwrSpec_Calc.(tmpname);

            file_internalmatch.Cd_CohPwr{sampcnt} = tmpdata.Cd_data; % Add to file_match struct
            file_internalmatch.Cohmean_CohPwr{sampcnt} = tmpdata.Cohrxx_mean_data;
            file_internalmatch.Cohd_CohPwr{sampcnt} = tmpdata.Cohrxx_d_data;
        end

        % At the point in the code: resampled calculations are collected
        % into cell and structure
        for ii = 1:Npop
            for jj = ii:Npop
                Cd_tmp_CohPwrResamp = zeros(length(freq),length(daxis));
                Coh_mean_tmp_CohPwrResamp = zeros(length(freq),1);
                Coh_d_tmp_CohPwrResamp = zeros(length(freq),length(daxis));
                %       Npair_tmp_CohPwrResamp = zeros(length(daxis),1);

                for mm = 1:resampcount
                    %  if size(file_match.C_d_Current{1,mm}{ii,jj}) == [length(daxis), 1] % 20x1
                    if ~isempty(file_internalmatch.Cd_CohPwr{1,mm}{ii,jj}) == 1
                        Cd_tmp_CohPwrResamp = Cd_tmp_CohPwrResamp + file_internalmatch.Cd_CohPwr{1,mm}{ii,jj};
                    end
                    if ~isempty(file_internalmatch.Cohmean_CohPwr{1,mm}{ii,jj}) == 1
                        Coh_mean_tmp_CohPwrResamp = Coh_mean_tmp_CohPwrResamp + file_internalmatch.Cohmean_CohPwr{1,mm}{ii,jj};
                    end
                    if ~isempty(file_internalmatch.Cohd_CohPwr{1,mm}{ii,jj}) == 1
                        Coh_d_tmp_CohPwrResamp = Coh_d_tmp_CohPwrResamp + file_internalmatch.Cohd_CohPwr{1,mm}{ii,jj};
                    end

                end

                Cd_mean_Internal{ii,jj} = Cd_tmp_CohPwrResamp/resampcount;
                Cohmean_Internal{ii,jj} = Coh_mean_tmp_CohPwrResamp/resampcount;
                Cohd_Internal{ii,jj} = Coh_d_tmp_CohPwrResamp/resampcount;
            end

        end

        ResampAvg.Cd_data = Cd_mean_Internal;
        ResampAvg.Cohrxx_mean_data = Cohmean_Internal;
        ResampAvg.Cohrxx_d_data = Cohd_Internal;
     %   ResampAvg.daxis = daxis;
     %   ResampAvg.Freq = freq;
      %  ResampAvg.eligibleToSamp = 

    elseif resampcount == 1
        sampcnt = resampcount;
        tmpname = sprintf('CohCalc_%s',num2str(sampcnt));
        ResampAvg = file.CohPwrSpec_Calc.(tmpname);

        ResampAvg.daxis = daxis;
        ResampAvg.Freq = freq;    
        
        ResampAvg.Npair_data = file.CohPwrSpec_Calc.Npair_data;
        ResampAvg.fs = file.CohPwrSpec_Calc.fs;
        ResampAvg.Twoverlap = file.CohPwrSpec_Calc.fs;
        ResampAvg.Tw = file.CohPwrSpec_Calc.Tw;
        ResampAvg.Nc = file.CohPwrSpec_Calc.Nc;

    end


    if contains(original_files(ff).name, 'ID') == 1
        %IndividData.CohPwrSpec_Calc{count} = ResampAvg;
        % IndividData.parmval{count} = file.ParamChange{param_rowm, param_coln};
           IndividData.CohPwrSpec_Calc{count} = file.CohPwrSpec_Calc;

        % !!!!!!!!!!!!!!!!!!!!!!!!!!
        IndividData.parmval{count} = file.paramVal1;% file.ParamChange{param_rowm, param_coln};
        %IndividData.parmval{count} = file.ParamChange{param_rowm, param_coln};
        IndividData.Npair{count} = file.CohPwrSpec_Calc.Npair_data;
        % IndividData.paramval(count) = file.ParamChange{param_rowm, param_coln};
        %CorrCurrentData_eachtrial.InputCurrent_Calc{count} = file.InputCurrent_Calc;

        count = count + 1;
    end
    %%%%%%% coherence data %%%%%
    cohpwspdata.Cd{ff} = ResampAvg.Cd_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x daxis
    cohpwspdata.Coh_mean{ff} = ResampAvg.Cohrxx_mean_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x 1
    cohpwspdata.Coh_d{ff} = ResampAvg.Cohrxx_d_data; %% upper tri cell Npop x Npop,  each input is length(f_range) x daxis
    % cohpwspdata.Freq = file.CohPwrSpec_Calc.Freq; %% frequency 0 - 500
    % cohpwspdata.daxis = file.CohPwrSpec_Calc.daxis; % vec: 1x20
    cohpwspdata.Npair{ff} = ResampAvg.Npair_data; %% upper tri cell Npop x Npop, 1 x daxis
    cohpwspdata.Nc{ff} = ResampAvg.Nc;


    %%%% NOTE: This needs to be changed to the correct parameter that
    %%%% changes across trials

     % !!!!!!!!!!!!!!!!!!!!!!!!!!
    data.alpha(ff) = file.paramVal1;%file.ParamChange{param_rowm, param_coln};
    %    data.alpha(ff) = file.ParamChange{param_rowm, param_coln};

end
save(matlabdata_filename_avgs, 'IndividData','-v7.3');
%save(matlabdata_filename_avgs, 'param_eachtrial','ParamChange_eachtrial','CorrCurrentData_eachtrial','CorrFRData_eachttrial','CohPwrSpec_Calc','-v7.3');
% save(matlabdata_filename_avgs, 'param_eachtrial','ParamChange_eachtrial','CorrCurrentData_eachtrial','CorrFRData_eachttrial','-v7.3');

%%%% save multicalcs
%Npop = 4;
daxis = file.CohPwrSpec_Calc.daxis;
freq = file.CohPwrSpec_Calc.Freq;
unique_param = unique(data.alpha); %finds unique values of changing parameter

Cd_data_trialavg = cell(1,length(unique_param));
Cohmean_data_trialavg = cell(1,length(unique_param));
Cohd_data_trialavg = cell(1,length(unique_param));
Npair_data_trialavg = cell(1,length(unique_param));
paramVal = [];


for ss = 1:length(unique_param)
    Cd_mean_CohPwr = cell(Npop,Npop);
    Cohmean_CohPwr = cell(Npop,Npop);
    Cohd_CohPwr = cell(Npop,Npop);
    Npair_CohPwr = cell(Npop,Npop);


    matches = data.alpha == unique_param(ss); % Grab data.alpha locations that match specific alpha
    match = 0; % Index to see the matches (running count)
    for kk = 1:length(matches) % Run through all possible matches to look for non zero matches
        if matches(kk) == 1 % If non zero match, add to file_match struct
            match = match + 1; % Increase index
            file_match.Cd_CohPwr{match} = cohpwspdata.Cd{kk}; % Add to file_match struct
            file_match.Cohmean_CohPwr{match} = cohpwspdata.Coh_mean{kk};
            file_match.Cohd_CohPwr{match} = cohpwspdata.Coh_d{kk};
            file_match.Npair_CohPwr{match} = cohpwspdata.Npair{kk};

        end
    end
    % At the point in the code: all the files that correspond to a single changing
    % parameter value (i.e different trials across the same
    % condition) are collected.  Next, calculate the averages
    % across all trials for one condition
    for ii = 1:Npop %compute mean correlation and cov across trials daxis [0,20]
        for jj = ii:Npop
            Cd_tmp_CohPwr = zeros(length(freq),length(daxis));
            Coh_mean_tmp_CohPwr = zeros(length(freq),1);
            Coh_d_tmp_CohPwr = zeros(length(freq),length(daxis));
            Npair_tmp_CohPwr = zeros(length(daxis),1);

            for mm = 1:match
                %  if size(file_match.C_d_Current{1,mm}{ii,jj}) == [length(daxis), 1] % 20x1
                if ~isempty(file_match.Cd_CohPwr{1,mm}{ii,jj}) == 1
                    Cd_tmp_CohPwr = Cd_tmp_CohPwr + file_match.Cd_CohPwr{1,mm}{ii,jj};
                end
                if ~isempty(file_match.Cohmean_CohPwr{1,mm}{ii,jj}) == 1
                    Coh_mean_tmp_CohPwr = Coh_mean_tmp_CohPwr + file_match.Cohmean_CohPwr{1,mm}{ii,jj};
                end
                if ~isempty(file_match.Cohd_CohPwr{1,mm}{ii,jj}) == 1
                    Coh_d_tmp_CohPwr = Coh_d_tmp_CohPwr + file_match.Cohd_CohPwr{1,mm}{ii,jj};
                end
                if ~isempty(file_match.Npair_CohPwr{1,mm}{ii,jj}) == 1
                    %  ss
                    %  ii
                    %  jj
                    %  file_match.Npair_CohPwr{1,mm}{ii,jj}
                    Npair_tmp_CohPwr = Npair_tmp_CohPwr + file_match.Npair_CohPwr{1,mm}{ii,jj};
                end

                %  end
            end
            Cd_mean_CohPwr{ii,jj} = Cd_tmp_CohPwr/match;
            Cohmean_CohPwr{ii,jj} = Coh_mean_tmp_CohPwr/match;
            Cohd_CohPwr{ii,jj} = Coh_d_tmp_CohPwr/match;
            Npair_CohPwr{ii,jj} = Npair_tmp_CohPwr/match;

        end
    end

    Cd_data_trialavg{1,ss} = Cd_mean_CohPwr;
    Cohmean_data_trialavg{1,ss} = Cohmean_CohPwr;
    Cohd_data_trialavg{1,ss} = Cohd_CohPwr;
    Npair_data_trialavg{1,ss} = Npair_CohPwr;
    paramVal(ss) = unique_param(ss);



end

TrialAvgData.Cd_data_trialavg = Cd_data_trialavg;
TrialAvgData.Cohmean_data_trialavg = Cohmean_data_trialavg;
TrialAvgData.Cohd_data_trialavg = Cohd_data_trialavg;
TrialAvgData.Npair_data_trialavg = Npair_data_trialavg;
TrialAvgData.param = paramVal;

save(matlabdata_filename_avgs, 'TrialAvgData','-append');






end



