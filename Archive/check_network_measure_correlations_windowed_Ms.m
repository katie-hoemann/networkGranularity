%% Instructions:
% run the following scripts:
% - assign_subject_communities_loop (set subjectFiles to 0)
% - check_subject_communities_loop (set subjectFiles to 0)
% - check_subject_clustering_loop (set subjectFiles to 0)
% - check_subject_participation_loop (set subjectFiles to 0)
% - check_subject_density_loop
% - check_subject_spatial_metrics_loop
% - granularity_calculations
% make sure the necessary .mat files are on your path

clear;
clc;

%% specify dataset
dataSet = 18; % set to 18, 39, or 88 and the rest will take care of itself!
print = 0; % set to 1 to generate spreadsheet with results

%% load data
filenames{1} = ['dataSet' num2str(dataSet) '_number_communities_windows.mat'];
%filenames{2} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_clustering.mat'];
%filenames{5} = ['dataSet' num2str(dataSet) '_participation.mat'];
%filenames{6} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_zICC.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end

%% specify subjects
if dataSet == 39 % remove subject 31 (row 29) for all calculations
    numCommunities(29,:) = [];
    %modularity(29,:) = [];
    modularity_VA(29,:) = [];
    clusterCoefGlobal(29,:) = [];
    %participation(29,:) = [];
    %gateway(29,:) = [];
    diversity(29,:) = [];
    networkDensity(29,:) = [];
    commRadius(29,:) = [];
    rtoZ(29,:) = [];
end
  
%% construct correlation matrix
% metrics = [numCommunities(:,2) modularity(:,2) modularity_VA(:,2) clusterCoefGlobal(:,2) participation(:,2) participation(:,5)...
%     gateway(:,2) gateway(:,5) diversity(:,2) diversity(:,5) networkDensity(:,5)];
metrics = [rtoZ(:,3) clusterCoefGlobal(:,2) networkDensity(:,5) numCommunities(:,2) diversity(:,2) diversity(:,5) commRadius(:,2) modularity_VA(:,2)];
idx = find(isnan(metrics));
[row,~] = ind2sub(size(metrics),idx);
metrics(row,:) = [];
metricMatrix = corrcoef(metrics);

%% organize and save output
correlations = array2table(metricMatrix);
% checks = {'numCommunities' 'modularity' 'modularityVA' 'clusteringCoef' 'posParticipation' 'negParticipation'...
%     'posGatewayNode' 'negGatewayNode' 'posDiversity' 'negDiversity' 'networkDensityW'}';
checks = {'zICC' 'clust' 'dens' 'numCom' 'posDiv' 'negDiv' 'comRad' 'modVB'}';
correlations = horzcat(checks,correlations);
% correlations.Properties.VariableNames = {'metric' 'numCommunities' 'modularity' 'modularityVA' 'clusteringCoef' 'posParticipation' 'negParticipation'...
%     'posGatewayNode' 'negGatewayNode' 'posDiversity' 'negDiversity' 'networkDensityW'};
correlations.Properties.VariableNames = {'measure' 'zICC' 'clusteringCoef' 'density' 'numCommunities' 'posDiversity' 'negDiversity' 'comRadius' 'modularityVB'};

if print == 1
    filename = ['dataSet' num2str(dataSet) '_polychoric_network_metric_correlations.xlsx'];
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
filename = ['dataSet' num2str(dataSet) '_heatmap_Pearson_measure_correlation_matrix'];
saveas(heatmap,filename,'tiff');

%% create scatter plots
figure;
metricsTable = array2table(metrics);
% metricsTable.Properties.VariableNames = {'numCom' 'mod' 'modVA' 'clust' 'posPa' 'negPa' 'posGa' 'negGa' 'posDi' 'negDi' 'dens'};
metricsTable.Properties.VariableNames = {'zICC' 'clust' 'dens' 'numCom' 'pDiv' 'nDiv' 'comRad' 'modVA'};
scatterplots = corrplot(metricsTable);
%scatterplots = plotmatrix(metrics); % another option    
filename = ['dataSet' num2str(dataSet) '_Pearson_measure_scatter_plots'];
% saveas(scatterplots,filename,'tiff');