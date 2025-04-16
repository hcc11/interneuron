
function Calc_CohPeakandThreshold(name,data_folder)
%%%%% This requires:
%%%% (1) Running calculation files: Calc_corrFR.m and Calc_cspd.m
%%%% (2) Generating summary files: DataExtraction_CohPwrPsc.m and
%%%% DataExtraction_FRCorr.m
%%%% saves 'CohPeakData' in ['AllCohPwr_Avgs_resamp_',name] 
Npop = 4;


tmpCohname = sprintf('AllCohPwr_Avgs_resamp_%s*',name);
tmpCohfile = dir([data_folder tmpCohname]);
tmpCohfile.name, 
file = load([data_folder tmpCohfile.name]);
matlabdata_filename = sprintf("%s%s",data_folder,tmpCohfile.name);


tmpFRname = sprintf('AllFRCorr_Avgs_%s*',name);
tmpFRfile = dir([data_folder tmpFRname]);
tmpFRfile.name
frfile = load([data_folder tmpFRfile.name]);

if isfield(file.TrialAvgData,'freq')==0 
    file.TrialAvgData.freq = 0:500;
end 

paramVal = sort(file.TrialAvgData.param);
for ii = 1:length(paramVal)
    %%%% get means of each pop and data from each file
    [x,y] = size(file.TrialAvgData.Cd_data_trialavg{1,ii});
    Cd_data(:,ii) = file.TrialAvgData.Cd_data_trialavg{1,ii}([1:x+1:y*x])';
    Cohmean_data(:,ii) = file.TrialAvgData.Cohmean_data_trialavg{1,ii}([1:x+1:y*x])';
    Cohd_data(:,ii) = file.TrialAvgData.Cohd_data_trialavg{1,ii}([1:x+1:y*x])';
    nuSum_data(ii,:) = frfile.CorrFRData_AcrossTrials.nuSim_data_trialavg{1,ii};
    Npair(ii,:) = file.TrialAvgData.Npair_data_trialavg{1,ii}([1:x+1:y*x])';
end
Cd_data = fixnan(Cd_data);
firstpeakcoh = zeros(Npop,length(paramVal));
firstpeakfreq = zeros(Npop, length(paramVal));
for pop = 1:Npop
    %pop
    for ii = 1:length(paramVal)
        %ii
        %tmpmeanCd = mean(Cd_data{pop,ii}(2:end,:),2);
        if sum(Npair{ii,pop})==0
            m = 0; I =1; 
        else
        tmpmeanCd = Cd_data{pop,ii}*Npair{ii,pop}/sum(Npair{ii,pop});
        [m, I] = max(tmpmeanCd);
        end
        firstpeakcoh(pop,ii) = m;
        firstpeakfreq(pop,ii) = file.TrialAvgData.freq(I);
        % if nnz(tmpmeanCd) ==
%         [cohpeak freqpeak] = findpeaks(tmpmeanCd);
%         if ~isempty(freqpeak) == 1
%             [~, indx1] = max(cohpeak(:,1));
%             firstpeakcoh(pop,ii) = cohpeak(indx1,:);
%             firstpeakfreq(pop,ii) = freqpeak(indx1,:);
%         else

%         end

    end
end
CohPeakData.firstPeak = firstpeakcoh;
CohPeakData.firstPeakFreq = firstpeakfreq;
CohPeakData.nuSim = nuSum_data;

save(matlabdata_filename,'CohPeakData','-append');

end