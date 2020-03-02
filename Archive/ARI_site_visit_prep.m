clear all;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
windowed = 0; % set to 1 to use data from sliding windows
SDs = 0; % set to 1 to use standard deviations from windowed data

%% load data
if windowed == 1
    filenames{1} = ['dataSet' num2str(dataSet) '_networkMeasureMeans_windows.mat'];
    filenames{2} = ['dataSet' num2str(dataSet) '_networkMeasureSDs_windows.mat'];
    filenames{3} = ['dataSet' num2str(dataSet) '_zICC.mat'];
    filenames{4} = ['dataSet' num2str(dataSet) '_affect.mat'];
    for i_file = 1:numel(filenames)
        load(filenames{i_file});
    end
else
    filenames{1} = ['dataSet' num2str(dataSet) '_communities.mat'];
    filenames{2} = ['dataSet' num2str(dataSet) '_modularity.mat'];
    filenames{3} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
    filenames{4} = ['dataSet' num2str(dataSet) '_clustering.mat'];
    filenames{5} = ['dataSet' num2str(dataSet) '_participation.mat'];
    filenames{6} = ['dataSet' num2str(dataSet) '_gateway.mat'];
    filenames{7} = ['dataSet' num2str(dataSet) '_diversity.mat'];
    filenames{8} = ['dataSet' num2str(dataSet) '_density.mat'];
    filenames{9} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
    filenames{10} = ['dataSet' num2str(dataSet) '_zICC.mat'];
    filenames{11} = ['dataSet' num2str(dataSet) '_affect.mat'];
    for i_file = 1:numel(filenames)
        load(filenames{i_file});
    end
end

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
TAStotal = rawData.data(:,17); %39 for session 2
RDEES = rawData.data(:,20)./14; %45 for session 2 % divide by 14 to get mean per item

% network metrics of interest
if windowed == 1
    if SDs == 1
        numCom = networkMeasureSDs(:,4);
        modVA = networkMeasureSDs(:,15);
        cluster = networkMeasureSDs(:,1);
        density = networkMeasureSDs(:,3);
        posDiv = networkMeasureSDs(:,10);
        negDiv = networkMeasureSDs(:,11);
    else
        numCom = networkMeasureMeans(:,4);
        modVA = networkMeasureMeans(:,15);
        cluster = networkMeasureMeans(:,1);
        density = networkMeasureMeans(:,3);
        posDiv = networkMeasureMeans(:,10);
        negDiv = networkMeasureMeans(:,11);
    end
else
    numCom = numCommunities(:,2); % number of communities at no threshold
    modVA = modularity_VA(:,2); % modularity for valence-assigned communities at no threshold
    cluster = clusterCoefGlobal(:,2); % global clustering at no threshold
    density = networkDensity(:,5); % weighted density at no threshold
    posDiv = diversity(:,2); % positive diversity at no threshold
    negDiv = diversity(:,5); % negative diversity at no threshold
end
    
%% other variables to test
zICC = rtoZ(:,3); % zICC total
mNeg = affect(:,2); % mean negative affect
sdNeg = affect(:,4); % std negative affect
mPos = affect(:,1); % mean positive affect
sdPos = affect(:,3); % std positive affect

%% construct correlation matrix
metrics = [numCom modVA cluster density posDiv negDiv];
metricNames = {'numCom' 'modVA' 'cluster' 'density' 'posDiv' 'negDiv'};
idx = find(isnan(metrics));
[row,~] = ind2sub(size(metrics),idx);
metrics(row,:) = [];
metricMatrix = corrcoef(metrics);

%% create heatmap
gran = zICC*-1;
measures = [gran metrics];
measureMatrix = corrcoef(measures);
labels = {'granularity' 'communities' 'valencePartition' 'clustering' 'density' 'posDiversity' 'negDiversity'}';

figure;
colormap('parula')
heatmap = imagesc(measureMatrix);
colorbar;
caxis([-1 1]);
set(gca,'XTick',1:1:11,'XTickLabel',labels);
xtickangle(45)
set(gca,'YTick',1:1:11,'YTickLabel',labels);
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_heatmap_measure_correlation_matrix_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_heatmap_measure_correlation_matrix_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_heatmap_measure_correlation_matrix'];
end
saveas(heatmap,filename,'tiff');

%% create scatter plots of significant relationships
figure;
scatter1 = scatter(numCom,mPHQ8total,[],[0 0.4470 0.7410],'filled');
if windowed == 1
    if SDs == 1
        xlabel('number of communities (SD)');
    else
        xlabel('number of communities (M)');
        xlim([1 5]);
        xticks([1 2 3 4 5]);
    end
else
    xlabel('number of communities');
    xlim([1 5]);
    xticks([1 2 3 4 5]);
end
ylim([0 24]);
yticks([4 8 12 16 20 24]);
ylabel('mean self-reported depression (PHQ8)');
h1 = lsline;
h1.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_numCom_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_numCom_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mPHQ8total_numCom_scatter_plot'];
end
saveas(scatter1,filename,'tiff');

figure;
scatter2 = scatter(modVA,mPHQ8total,[],[0.8500 0.3250 0.0980],'filled');
if windowed == 1
    if SDs == 1
        xlabel('fit of valence-based partition (SD)');
    else
        xlabel('fit of valence-based partition (M)');
        xlim([0.2 0.8]);
    end
else
    xlabel('fit of valence-based partition');
    xlim([0.2 0.8]);
end
ylim([0 24]);
yticks([4 8 12 16 20 24]);
ylabel('mean self-reported depression (PHQ8)');
h2 = lsline;
h2.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_modVA_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_modVA_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mPHQ8total_modVA_scatter_plot'];
end
saveas(scatter2,filename,'tiff');

figure;
scatter3 = scatter(posDiv,mPHQ8total,[],[0.4660 0.6740 0.1880],'filled');
if windowed == 1
    if SDs == 1
        xlabel('diversity of positive connections between communities (SD)');
        xlim([0 0.3]);
    else
        xlabel('diversity of positive connections between communities (M)');
        xlim([0 0.7]);
    end
else
    xlabel('diversity of positive connections between communities');
    xlim([0 0.7]);
end
ylim([0 24]);
yticks([4 8 12 16 20 24]);
ylabel('mean self-reported depression (PHQ8)');
h3 = lsline;
h3.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_posDiv_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mPHQ8total_posDiv_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mPHQ8total_posDiv_scatter_plot'];
end
saveas(scatter3,filename,'tiff');

figure;
scatter4 = scatter(numCom,mGADtotal,[],[0 0.4470 0.7410],'filled');
if windowed == 1
    if SDs == 1
        xlabel('number of communities (SD)');
    else
        xlabel('number of communities (M)');
        xlim([1 5]);
        xticks([1 2 3 4 5]);
    end
else
    xlabel('number of communities');
    xlim([1 5]);
    xticks([1 2 3 4 5]);
end
ylim([0 21]);
yticks([3 6 9 12 15 18 21]);
ylabel('mean self-reported anxiety(GAD)');
h1 = lsline;
h1.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_numCom_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_numCom_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mGADtotal_numCom_scatter_plot'];
end
saveas(scatter4,filename,'tiff');

figure;
scatter5 = scatter(modVA,mGADtotal,[],[0.8500 0.3250 0.0980],'filled');
if windowed == 1
    if SDs == 1
        xlabel('fit of valence-based partition (SD)');
    else
        xlabel('fit of valence-based partition (M)');
        xlim([0.2 0.8]);
    end
else
    xlabel('fit of valence-based partition');
    xlim([0.2 0.8]);
end
ylim([0 21]);
yticks([3 6 9 12 15 18 21]);
ylabel('mean self-reported anxiety (GAD)');
h2 = lsline;
h2.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_modVA_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_modVA_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mGADtotal_modVA_scatter_plot'];
end
saveas(scatter5,filename,'tiff');

figure;
scatter6 = scatter(posDiv,mGADtotal,[],[0.4660 0.6740 0.1880],'filled');
if windowed == 1
    if SDs == 1
        xlabel('diversity of positive connections between communities (SD)');
        xlim([0 0.3]);
    else
        xlabel('diversity of positive connections between communities (M)');
        xlim([0 0.7]);
    end
else
    xlabel('diversity of positive connections between communities');
    xlim([0 0.7]);
end
ylim([0 21]);
yticks([3 6 9 12 15 18 21]);
ylabel('mean self-reported anxiety (GAD)');
h3 = lsline;
h3.Color = 'k';
if windowed == 1
    if SDs == 1
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_posDiv_scatter_plot_window_SDs'];
    else
        filename = ['dataSet' num2str(dataSet) '_mGADtotal_posDiv_scatter_plot_window_Ms'];
    end
else
    filename = ['dataSet' num2str(dataSet) '_mGADtotal_posDiv_scatter_plot'];
end
saveas(scatter6,filename,'tiff');
    