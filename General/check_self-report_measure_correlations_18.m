clear;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
windowed = 0; % set to 1 to use data from sliding windows
SDs = 0; % set to 1 to use standard deviations from windowed data
print = 0; % set to 1 to generate spreadsheet with results

%% load self-report data
rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASItotal = rawData.data(:,7)-21; % subtract 21 to bring in line with original rating scale
ERStotal = rawData.data(:,11); 
GADtotal_s1 = rawData.data(:,12); % 12 for session 1; 34 for session 2
GADtotal_s2 = rawData.data(:,34);
mGADtotal = (GADtotal_s1+GADtotal_s2)/2;
PSStotal = rawData.data(:,21).*3.5; %46 for session 2 % multiply by 3.5 to bring in line with norms from 14-item measure
PHQ15total = rawData.data(:,22); %47 for session 2
PHQ8total_s1 = rawData.data(:,23); %48 for session 2
PHQ8total_s2 = rawData.data(:,48);
mPHQ8total = (PHQ8total_s1+PHQ8total_s2)/2;
neuroticism = rawData.data(:,49);
TAStotal_s1 = rawData.data(:,17); %39 for session 2
TAStotal_s2 = rawData.data(:,39); 
mTAStotal = (TAStotal_s1+TAStotal_s2)/2;
RDEES = rawData.data(:,20)./14; %45 for session 2 % divide by 14 to get mean per item

%% load non-granularity derived measures
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_measureSDs_windows.mat'];
        load(filename);
        zICC = measureSDs(:,3); % zICC total
        instab = measureSDs(:,5); % emotional instability
        mPos = measureSDs(:,9); % mean positive affect
        mNeg = measureSDs(:,10); % mean negative affect
        sdPos = measureSDs(:,11); % std positive affect
        sdNeg = measureSDs(:,12); % std negative affect
    else
        filename = ['dataSet' num2str(dataSet) '_measureMeans_windows.mat'];
        load(filename);
        zICC = measureMeans(:,3); % zICC total
        instab = measureMeans(:,5); % emotional instability
        mPos = measureMeans(:,9); % mean positive affect
        mNeg = measureMeans(:,10); % mean negative affect
        sdPos = measureMeans(:,11); % std positive affect
        sdNeg = measureMeans(:,12); % std negative affect
    end
else
    filename = ['dataSet' num2str(dataSet) '_measures.mat'];
    load(filename);
    zICC = measures(:,3); % zICC total
    instab = measures(:,5); % emotional instability
    mPos = measures(:,9); % mean positive affect
    mNeg = measures(:,10); % mean negative affect
    sdPos = measures(:,11); % std positive affect
    sdNeg = measures(:,12); % std negative affect
end
  
%% construct correlation matrix of self-report measures
measures = [ASItotal(:) ERStotal(:) GADtotal_s1(:) PSStotal(:) PHQ15total(:) PHQ8total_s1(:) neuroticism(:) TAStotal_s1(:) RDEES(:)];
idx = find(isnan(measures));
[row,~] = ind2sub(size(measures),idx);
measures(row,:) = [];
measureMatrix = corrcoef(measures);

%% organize and save output
correlations = array2table(measureMatrix);
questionnaires = {'ASI3' 'ERS' 'GAD7' 'PSS' 'PHQ15' 'PHQ8' 'neuroticism' 'TAS20' 'RDEES'}';
correlations = horzcat(questionnaires,correlations);
correlations.Properties.VariableNames = {'measure' 'ASI3' 'ERS' 'GAD7' 'PSS' 'PHQ15' 'PHQ8' 'neuroticism' 'TAS20' 'RDEES'}';

if print == 1
    filename = ['dataSet' num2str(dataSet) '_self-report_measure_correlations.xlsx'];
    writetable(correlations,filename);
end

%% create heatmap
figure;
colormap('parula')
heatmap = imagesc(measureMatrix);
colorbar;
caxis([-1 1]);
set(gca,'XTick',1:1:11,'XTickLabel',questionnaires);
xtickangle(45)
set(gca,'YTick',1:1:11,'YTickLabel',questionnaires);
filename = ['dataSet' num2str(dataSet) '_heatmap_self-report_measure_correlation_matrix'];
saveas(heatmap,filename,'tiff');

%% create histograms
figure;
histogram1 = histogram(ASItotal(:));
xlim([0 72]);
ylim([0 50]);
vline([13 29],{'k:','k:'},{'Non-clinical M','Clinical M'});
title('ASI-3 scores');
filename = ['dataSet' num2str(dataSet) '_ASI3_scores'];
saveas(histogram1,filename,'tiff');

figure;
histogram2 = histogram(ERStotal(:));
xlim([0 84]);
ylim([0 50]);
vline([30 45],{'k:','k:'},{'Non-clinical M','Clinical M'});
title('ERS scores');
filename = ['dataSet' num2str(dataSet) '_ERS_scores'];
saveas(histogram2,filename,'tiff');

figure;
histogram3 = histogram(mGADtotal(:));
xlim([0 21]);
ylim([0 50]);
vline([5 10 15],{'g','b','r'},{'mild','moderate','severe'});
title('Mean GAD-7 scores');
filename = ['dataSet' num2str(dataSet) '_mGAD7_scores'];
saveas(histogram3,filename,'tiff');

figure;
histogram4 = histogram(PSStotal(:));
xlim([0 56]);
ylim([0 50]);
vline([0 14 27],{'g','b','r'},{'low','moderate','high'});
title('PSS scores');
filename = ['dataSet' num2str(dataSet) '_PSS_scores'];
saveas(histogram4,filename,'tiff');

figure;
histogram5 = histogram(PHQ15total(:));
xlim([0 30]);
ylim([0 50]);
vline([5 10 15],{'g','b','r'},{'low','medium','high'});
title('PHQ-15 scores');
filename = ['dataSet' num2str(dataSet) '_PHQ15_scores'];
saveas(histogram5,filename,'tiff');

figure;
histogram6 = histogram(mPHQ8total(:));
xlim([0 24]);
ylim([0 50]);
vline([5 10 15 20],{'g','b','k','r'},{'mild','moderate','mod severe','severe'});
title('Mean PHQ-8 scores');
filename = ['dataSet' num2str(dataSet) '_mPHQ8_scores'];
saveas(histogram6,filename,'tiff');

figure;
histogram7 = histogram(neuroticism(:));
xlim([20 180]);
ylim([0 50]);
vline([55 75 95 115],{'g','b','k','r'},{'low','average','high','very high'});
title('Neuroticism scores');
filename = ['dataSet' num2str(dataSet) '_neuroticism_scores'];
saveas(histogram7,filename,'tiff');

figure;
histogram8 = histogram(mTAStotal(:));
xlim([20 100]);
ylim([0 50]);
vline([20 51 61],{'g','b','r'},{'none','possible','alexithymia'});
title('Mean TAS-20 scores');
filename = ['dataSet' num2str(dataSet) '_mTAS20_scores'];
saveas(histogram8,filename,'tiff');

figure;
histogram9 = histogram(RDEES(:));
xlim([1 5]);
ylim([0 50]);
vline(4,'k:','Non-clinical M');
title('RDEES scores')
filename = ['dataSet' num2str(dataSet) '_RDEES_scores'];
saveas(histogram9,filename,'tiff');

figure;
histogram10 = histogram(instab(:));
ylim([0 50]);
title('Emotional instability')
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_instability_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_instability_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_instability'];
end
saveas(histogram10,filename,'tiff');

figure;
histogram11 = histogram(mNeg(:));
ylim([0 50]);
title('Mean negative affect')
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mNeg_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mNeg_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mNeg'];
end
saveas(histogram11,filename,'tiff');

figure;
histogram12 = histogram(sdNeg(:));
ylim([0 50]);
title('Standard deviation negative affect')
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_sdNeg_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_sdNeg_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_sdNeg'];
end
saveas(histogram12,filename,'tiff');

figure;
histogram13 = histogram(zICC(:));
ylim([0 50]);
title('zICC (granularity)')
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_zICC_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_zICC_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_zICC'];
end
saveas(histogram13,filename,'tiff');

%% create scatter plots
figure;
measuresTable = array2table(measures);
measuresTable.Properties.VariableNames = {'ASI3' 'ERS' 'GAD7' 'PSS' 'PHQ15' 'PHQ8' 'neuroticism' 'TAS20' 'RDEES'};
scatterplots = corrplot(measuresTable);
% scatterplots = plotmatrix(measures); % another option    
% filename = ['dataSet' num2str(dataSet) '_self-report_measures_scatter_plots'];
% saveas(scatterplots,filename,'tiff'); % not working

close all;