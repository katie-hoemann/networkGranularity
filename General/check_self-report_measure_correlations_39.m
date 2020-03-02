clear;
clc;

%% specify dataset
dataSet = 39; % must be set to 39
windowed = 1; % set to 1 to use data from sliding windows
SDs = 1; % set to 1 to use standard deviations from windowed data
print = 0; % set to 1 to generate spreadsheet with results

%% load criterion validity measures
rawData = importdata('39ESdata_QuestionnaireData.xlsx');
BAItotal = rawData.data(:,9); 
BDItotal = rawData.data(:,12);
SWLStotal = rawData.data(:,13);
TAStotal = rawData.data(:,2);
RDEES = rawData.data(:,6);

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

%% remove subject 31 (row 29) for all calculations
BAItotal(29,:) = [];
BDItotal(29,:) = [];
SWLStotal(29,:) = [];
TAStotal(29,:) = [];
RDEES(29,:) = [];
zICC(29,:) = [];
instab(29,:) = [];
mPos(29,:) = [];
mNeg(29,:) = [];
sdPos(29,:) = [];
sdNeg(29,:) = [];
   
%% construct correlation matrix
measures = [BAItotal(:) BDItotal(:) SWLStotal(:) TAStotal(:) RDEES(:)];
idx = find(isnan(measures));
[row,~] = ind2sub(size(measures),idx);
measures(row,:) = [];
measureMatrix = corrcoef(measures);

%% organize and save output
correlations = array2table(measureMatrix);
questionnaires = {'BAI' 'BDI' 'SWLS' 'TAS20' 'RDEES'}';
correlations = horzcat(questionnaires,correlations);
correlations.Properties.VariableNames = {'measure' 'BAI' 'BDI' 'SWLS' 'TAS20' 'RDEES'}';

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
histogram1 = histogram(BAItotal(:));
xlim([0 63]);
ylim([0 76]);
vline([0 10 19 30],{'g','b','k','r'},{'normal','mild','moderate','severe'});
title('BAI scores');
filename = ['dataSet' num2str(dataSet) '_BAI_scores'];
saveas(histogram1,filename,'tiff');

figure;
histogram2 = histogram(BDItotal(:));
xlim([0 63]);
ylim([0 76]);
vline([0 10 19 30],{'g','b','k','r'},{'normal','mild','moderate','severe'});
title('BDI scores');
filename = ['dataSet' num2str(dataSet) '_BDI_scores'];
saveas(histogram2,filename,'tiff');

figure;
histogram3 = histogram(SWLStotal(:));
xlim([5 35]);
ylim([0 76]);
vline(23.5,'k:','Non-clinical M');
title('SWLS scores');
filename = ['dataSet' num2str(dataSet) '_SWLS_scores'];
saveas(histogram3,filename,'tiff');

figure;
histogram4 = histogram(TAStotal(:));
xlim([20 100]);
ylim([0 76]);
vline([20 51 61],{'g','b','r'},{'none','possible','alexithymia'});
title('TAS-20 scores');
filename = ['dataSet' num2str(dataSet) '_TAS20_scores'];
saveas(histogram4,filename,'tiff');

figure;
histogram5 = histogram(RDEES(:));
xlim([0 20]);
ylim([0 76]);
vline(4,'k:','Non-clinical M');
title('RDEES scores');
filename = ['dataSet' num2str(dataSet) '_RDEES_scores'];
saveas(histogram5,filename,'tiff');

figure;
histogram6 = histogram(instab(:));
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
saveas(histogram6,filename,'tiff');

figure;
histogram7 = histogram(mNeg(:));
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
saveas(histogram7,filename,'tiff');

figure;
histogram8 = histogram(sdNeg(:));
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
saveas(histogram8,filename,'tiff');

figure;
histogram9 = histogram(zICC(:));
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
saveas(histogram9,filename,'tiff');

%% create scatter plots
figure;
measureTable = array2table(measures);
measureTable.Properties.VariableNames = {'BAI' 'BDI' 'SWLS' 'TAS20' 'RDEES'};
scatterplots = corrplot(measureTable);
% scatterplots = plotmatrix(measures); % another option    
% filename = ['dataSet' num2str(dataSet) '_self-report_measures_scatter_plots'];
% saveas(scatterplots,filename,'tiff'); % not working 

close all;