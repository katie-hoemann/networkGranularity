%% Instructions:
% make sure the necessary .mat files are on your path

clear;
clc;

%% specify dataset
dataSet = 39; % set to 18, 39, or 88 and the rest will take care of itself!
windowed = 1; % set to 1 to use data from sliding windows
SDs = 0; % set to 1 to use standard deviations from windowed data
print = 0; % set to 1 to generate spreadsheet with results

%% load data
if windowed == 1
    if SDs == 1     
        filenames{1} = ['dataSet' num2str(dataSet) '_networkMeasureSDs_Pearson_windows.mat'];
        filenames{2} = ['dataSet' num2str(dataSet) '_zICC.mat'];
        for i_file = 1:numel(filenames)
            load(filenames{i_file});
        end
    else
        filenames{1} = ['dataSet' num2str(dataSet) '_networkMeasureMeans_Pearson_windows.mat'];
        filenames{2} = ['dataSet' num2str(dataSet) '_zICC.mat'];
        for i_file = 1:numel(filenames)
            load(filenames{i_file});
        end
    end
else
    filenames{1} = ['dataSet' num2str(dataSet) '_networkMeasures_Pearson.mat'];
    filenames{2} = ['dataSet' num2str(dataSet) '_zICC.mat'];
    for i_file = 1:numel(filenames)
        load(filenames{i_file});
    end
end

%% specify subjects
if dataSet == 39 % remove subject 31 (row 29) for all calculations
    if windowed == 1
        if SDs == 1
            networkMeasureSDs(29,:) = [];
        else
            networkMeasureMeans(29,:) = [];
        end
    else
        networkMeasures(29,:) = [];
    end
    rtoZ(29,:) = [];
end

%% specify network metrics of interest
if windowed == 1
    if SDs == 1
        cluster = networkMeasureSDs(:,1);
        density = networkMeasureSDs(:,3);
        numCom = networkMeasureSDs(:,4);
        posDiv = networkMeasureSDs(:,10);
        negDiv = networkMeasureSDs(:,11);
        comRad = networkMeasureSDs(:,14);
        modVB = networkMeasureSDs(:,15);
    else
        cluster = networkMeasureMeans(:,1);
        density = networkMeasureMeans(:,3);
        numCom = networkMeasureMeans(:,4);
        posDiv = networkMeasureMeans(:,10);
        negDiv = networkMeasureMeans(:,11);
        comRad = networkMeasureMeans(:,14);
        modVB = networkMeasureMeans(:,15);
    end
else
    cluster = networkMeasures(:,1);
    density = networkMeasures(:,3);
    numCom = networkMeasures(:,4);
    posDiv = networkMeasures(:,10);
    negDiv = networkMeasures(:,11);
    comRad = networkMeasures(:,14);
    modVB = networkMeasures(:,15);
end
zICC = rtoZ(:,3);
  
%% construct correlation matrix
metrics = [zICC cluster density numCom posDiv negDiv comRad modVB];
idx = find(isnan(metrics));
[row,~] = ind2sub(size(metrics),idx);
metrics(row,:) = [];
metricMatrix = corrcoef(metrics);

%% organize and save output
correlations = array2table(metricMatrix);
checks = {'zICC' 'clust' 'dens' 'nCom' 'posDiv' 'negDiv' 'comR' 'modVB'}';
correlations = horzcat(checks,correlations);
correlations.Properties.VariableNames = {'measure' 'zICC' 'clusteringCoef' 'density' 'numCommunities' 'posDiversity' 'negDiversity' 'comRadius' 'modularityVB'};

if print == 1
    if windowed == 1
        if SDs == 1
            filename = ['dataSet' num2str(dataSet) '_Pearson_measure_SDs_correlations.xlsx'];
        else
            filename = ['dataSet' num2str(dataSet) '_Pearson_measure_Ms_correlations.xlsx'];
        end
    else
        filename = ['dataSet' num2str(dataSet) '_Pearson_measure_correlations.xlsx'];
    end
    writetable(correlations,filename);
end

%% create heatmap
figure;
colormap('parula')
heatmap = imagesc(metricMatrix);
colorbar;
caxis([-1 1]);
set(gca,'XTick',1:1:11,'XTickLabel',checks);
xtickangle(45)
set(gca,'YTick',1:1:11,'YTickLabel',checks);
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_heatmap_Pearson_measure_SDs_correlation_matrix'];
    else
        filename = ['dataSet' num2str(dataSet) '_heatmap_Pearson_measure_Ms_correlation_matrix'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_heatmap_Pearson_measure_correlation_matrix'];
end
saveas(heatmap,filename,'tiff');

%% create scatter plots
figure;
metricsTable = array2table(metrics);
metricsTable.Properties.VariableNames = {'zICC' 'clust' 'dens' 'nCom' 'pDiv' 'nDiv' 'comR' 'modVA'};
[R,P,scatterplots] = corrplot(metricsTable);
title = findall(scatterplots(1).Parent.Parent, 'type', 'text', 'String', '{\bf Correlation Matrix}'); 
title.String = [];
%scatterplots = plotmatrix(metrics); % another option
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_Pearson_measure_SDs_scatter_plots'];
    else
        filename = ['dataSet' num2str(dataSet) '_Pearson_measure_Ms_scatter_plots'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_Pearson_measure_scatter_plots'];
end
%saveas(scatterplots,filename,'tiff'); % not working for some reason