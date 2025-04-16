%%%%%% pull data from trial runs and calculations into two matlab files:
%%%%%%%   
%%%%%%%         uses 'paramVal1' in each file as param values. 
function DataExtraction_FRCorr(data_folder, filename, like_filename)
%%%%%%% saves a file as ['AllFRCorr_Avgs_', filename] in data_folder 

%if ~isfield(option, 'dataextract') option.dataextract=0; end
%if option.dataextract == 1
    % % Initial setup: load files
    %like_newfilename = erase(like_filename, 'ID1')
    like_newfilename = like_filename;
    
    original_files = dir([data_folder like_newfilename])
    file_num = length(original_files)


   % matlabdata_filename_individ = sprintf("%s%s%s",data_folder,'AllData_Individ_', filename);
    matlabdata_filename_avgs = sprintf("%s%s%s",data_folder,'AllFRCorr_Avgs_', filename);
    %%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%
    %%%%%%%%%   %%%%%%%%%    %%%%%%%%%    %%%%%%%%%

    tic;
    %if option.acrosstrialcalcs == 1
    disp('FR Corr Data Grab')
    count = 1
    for ff = 1:file_num

        ff
        % % load .mat simulation file
        original_files(ff).name
        file = load([data_folder original_files(ff).name]);

        if contains(original_files(ff).name, 'ID') == 1
           IndividData.corrFR_Calc{count} = file.corrCalc;

           % %
           %%%% static input {param_rowm, param_coln} == 12, 2
          % IndividData.paramval(count) = file.ParamChange{param_rowm, param_coln};%file.param.Iapp(popwithinput,1);%file.ParamChange{param_rowm, param_coln};
           IndividData.paramval(count) = file.paramVal1;           
           %%%% if heterogenous input 
         %  IndividData.paramval(count) = file.ParamChange{13, 2}(popwithinput);

           
           IndividData.param = file.param;
           IndividData.nuSim{count} = file.nuSim;
           IndividData.s1{count} = file.s1;
           count = count + 1;
        end


        %%%% local storage for averaging across trials %%%%
        %%% store data for each file to be averaged by trial count

        %%%%%%%% correlation data on FR %%%%%%%
        corrdata.chi{ff} = file.corrCalc.chi;
        corrdata.freq{ff} = file.corrCalc.freq;
        corrdata.psdx{ff} = file.corrCalc.psdx;
        corrdata.psdx_peaks{ff} = file.corrCalc.peakdata;

        %%% store data for each file to be averaged by trial count
        corrdata.C_d{ff} = file.corrCalc.C_d; % 1x nfiles structure
        corrdata.COV_d{ff}=file.corrCalc.COV_d;
        corrdata.fano{ff} = file.corrCalc.Var./file.corrCalc.SC';
        corrdata.fanomean{ff} = file.corrCalc.fanomean;
        corrdata.Cbar{ff} = file.corrCalc.Cbar;
        corrdata.COVbar{ff} = file.corrCalc.COVbar;
        corrdata.Var{ff} = file.corrCalc.Var;
        corrdata.SC{ff} = file.corrCalc.SC;

        corrdata.nuSim{ff} = file.nuSim; %average firing rates recorded
        %%%% NOTE: This needs to be changed to the correct parameter that
        %%%% changes across trials

        %%%% static == {param_rowm, param_coln} = {12,2}
        %%%%heterogen === {param_rowm, param_coln} = {13,2}
       % data.alpha(ff) = file.ParamChange{param_rowm, param_coln};
       % data.alpha(ff) = file.ParamChange{13, param_coln}(popwithinput);

        % !!!!!!!!!!!!!!!
        data.alpha(ff) = file.paramVal1;

    end
    %save(matlabdata_filename_individ, 'CorrCurrentData_eachtrial','CorrFRData_eachttrial','-v7.3');
    save(matlabdata_filename_avgs, 'IndividData','-v7.3');
   % save(matlabdata_filename_avgs, 'param_eachtrial','ParamChange_eachtrial','CorrCurrentData_eachtrial','CorrFRData_eachttrial','-v7.3');

    %%%% save multicalcs
    Npop = 4;
    daxis=file.corrCalc.daxis;
    unique_param = unique(data.alpha); %finds unique values of changing paramete

    nuSim_mean_Corr = zeros(length(unique_param),Npop);
    chi_mean_Corr = zeros(length(unique_param),Npop);
    psdx_mean_Corr = cell(1,length(unique_param));
    psdxPeak_mean_Corr = zeros(length(unique_param),Npop);


    for ss = 1:length(unique_param)
        Cbar_mean_Corr = cell(Npop,Npop);
        C_d_mean_Corr = cell(Npop,Npop);
        COV_d_mean_Corr = cell(Npop,Npop);

        matches = data.alpha == unique_param(ss); % Grab data.alpha locations that match specific alpha
        match = 0; % Index to see the matches (running count)
        for kk = 1:length(matches) % Run through all possible matches to look for non zero matches
            if matches(kk) == 1 % If non zero match, add to file_match struct
                match = match + 1; % Increase index
                file_match.Cbar_Corr{match} = corrdata.Cbar{kk}; % Add to file_match struct
                file_match.C_d_Corr{match} = corrdata.C_d{kk};
                file_match.COV_d_Corr{match} = corrdata.COV_d{kk};
                file_match.nuSim_Corr{match} = corrdata.nuSim{kk};
                file_match.psdx_Corr{match} = corrdata.psdx{kk};
                file_match.psdx_peaks_Corr{match} = corrdata.psdx_peaks{kk};
                file_match.chi_Corr{match} = corrdata.chi{kk};
            end
        end
        % At the point in the code: all the files that correspond to a single changing
        % parameter value (i.e different trials across the same
        % condition) are collected.  Next, calculate the averages
        % across all trials for one condition
        for ii = 1:Npop %compute mean correlation and cov across trials daxis [0,20]
            for jj = ii:Npop
                Cbar_tmp_data_Corr = zeros(match,1);
                Ctmp_data_Corr = zeros(match,length(daxis));
                COVtmp_data_Corr = zeros(match,length(daxis));

                for mm = 1:match
                    Cbar_tmp_data_Corr(mm,1) = file_match.Cbar_Corr{1,mm}(ii,jj);
                    Ctmp_data_Corr(mm,:) = file_match.C_d_Corr{1,mm}{ii,jj};
                    COVtmp_data_Corr(mm,:) = file_match.COV_d_Corr{1,mm}{ii,jj};
                end
                
                Cbar_mean_Corr{ii,jj} = mean(Cbar_tmp_data_Corr);
                C_d_mean_Corr{ii,jj} = mean(Ctmp_data_Corr); % calculate mean across trials
                COV_d_mean_Corr{ii,jj} = mean(COVtmp_data_Corr);
            end
        end

        %%% collect firing rate averages in each file, average across
        %%% trials
        nuSimtmp_data_Corr = zeros(length(file_match.nuSim_Corr),Npop);
        for mm = 1:length(file_match.nuSim_Corr)
            for jj = 1:Npop
                nuSimtmp_data_Corr(mm,jj) = file_match.nuSim_Corr{1,mm}(1,jj);
            end
        end
        if size(nuSimtmp_data_Corr,1) == 1
            nuSim_mean_Corr(ss,:) = mean(nuSimtmp_data_Corr,1);
        else
            nuSim_mean_Corr(ss,:) = mean(nuSimtmp_data_Corr);
        end

        %%% collect powerspec (FFT) measurment averages in each file, average across
        %%% trials
        %  if option.powerspec_trialavg == 1
        pwrspectmp_data_Corr = zeros(length(file_match.psdx_Corr{1,1}),Npop);
        for mm = 1:length(file_match.psdx_Corr)
            pwrspectmp_data_Corr = pwrspectmp_data_Corr + file_match.psdx_Corr{1,mm};
        end

        if length(file_match.psdx_Corr) == 1
            pwrspec_mean_Corr{1,ss} = pwrspectmp_data_Corr;
        else
            pwrspec_mean_Corr{1,ss} = pwrspectmp_data_Corr/length(file_match.psdx_Corr);
        end

        % end
        peaktmp_data_Corr = zeros(length(file_match.psdx_peaks_Corr),Npop);
        for mm = 1:length(file_match.psdx_peaks_Corr)
            for jj = 1:Npop
                peaktmp_data_Corr(mm,jj) = file_match.psdx_peaks_Corr{1,mm}(1,jj);
            end
        end
        if size(peaktmp_data_Corr,1) == 1
            psdxPeak_mean_Corr(ss,:) = mean(peaktmp_data_Corr,1);
        else
            psdxPeak_mean_Corr(ss,:) = mean(peaktmp_data_Corr);
        end

        %%% collect synchony measurment averages in each file, average across
        %%% trials
        chitmp_data_Corr = zeros(length(file_match.chi_Corr),Npop);
        for mm = 1:length(file_match.chi_Corr)
            for jj = 1:Npop
                chitmp_data_Corr(mm,jj) = file_match.chi_Corr{1,mm}(1,jj);
            end
        end

        if size(chitmp_data_Corr,1) == 1
            chi_mean_Corr(ss,:) = mean(chitmp_data_Corr,1);
        else
            chi_mean_Corr(ss,:) = mean(chitmp_data_Corr);
        end

        %%% structure of plotdataset: row correspond to unique(unique_param); col
        %%% correspond to [C_d_mean, C_d_std, COV_d_mean, COV_d_std,...]
        C_d_mean_data_trialavg{1,ss} = C_d_mean_Corr;
        COV_d_data_trialavg{1,ss} = COV_d_mean_Corr;
        Cbar_data_trialavg{1,ss} = Cbar_mean_Corr;
        nuSim_data_trialavg{1,ss} = nuSim_mean_Corr(ss,:);
        pwrspec_data_trialavg{1,ss} = pwrspec_mean_Corr{1,ss};
        psdxPeak_data_trialavg{1,ss} = psdxPeak_mean_Corr(ss,:);
        chi_data_trialavg{1,ss} = chi_mean_Corr(ss,:);
        paramVal(ss)= unique_param(ss);

    end

CorrFRData_AcrossTrials.C_d_mean_data_trialavg = C_d_mean_data_trialavg;
CorrFRData_AcrossTrials.COV_d_data_trialavg = COV_d_data_trialavg;
CorrFRData_AcrossTrials.Cbar_data_trialavg = Cbar_data_trialavg;
CorrFRData_AcrossTrials.nuSim_data_trialavg = nuSim_data_trialavg;
CorrFRData_AcrossTrials.pwrspec_data_trialavg = pwrspec_data_trialavg;
CorrFRData_AcrossTrials.psdxPeak_data_trialavg = psdxPeak_data_trialavg;
CorrFRData_AcrossTrials.chi_data_trialavg = chi_data_trialavg;
CorrFRData_AcrossTrials.paramVal = paramVal;

save(matlabdata_filename_avgs, 'CorrFRData_AcrossTrials', '-append');


end



