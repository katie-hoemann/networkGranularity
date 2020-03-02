dataSet = 18;
load(['dataSet' num2str(dataSet) '_networkMeasures.mat']);
load(['dataSet' num2str(dataSet) '_networkMeasureMeans_windows.mat']);
load(['dataSet' num2str(dataSet) '_networkMeasureSimpleMeans_windows.mat']);
load(['dataSet' num2str(dataSet) '_networkMeasureWeightedMeans_windows.mat']);
for i_var = 1:size(networkMeasures,2)
    varCorr(1,i_var) = corr(networkMeasures(:,i_var),networkMeasureMeans(:,i_var));
    %varCorr(2,i_var) = corr(networkMeasures(:,i_var),networkMeasureSimpleMeans(:,i_var));
    %varCorr(3,i_var) = corr(networkMeasures(:,i_var),networkMeasureWeightedMeans(:,i_var));
    %varCorr(4,i_var) = corr(networkMeasureMeans(:,i_var),networkMeasureSimpleMeans(:,i_var));
    %varCorr(5,i_var) = corr(networkMeasureMeans(:,i_var),networkMeasureWeightedMeans(:,i_var));
    %varCorr(6,i_var) = corr(networkMeasureSimpleMeans(:,i_var),networkMeasureWeightedMeans(:,i_var));
end
variableNames = {'cluster','clusterSD','density','numCom','mod','pPos','pNeg','pPosSD','pNegSD','dPos','dNeg','dPosSD','dNegSD','comRadius','modVA'};
%rowNames = {'static v aggM','static v simpleM','static v weightM','aggM v simpleM','aggM v weightM','simpleM v weightM'}';
rowNames = {'static v aggM'};
varCorr_Table = array2table(varCorr,'VariableNames',variableNames);
varCorr_Table = horzcat(rowNames,varCorr_Table);
varCorr_Table.Properties.VariableNames{1} = 'comparison';
writetable(varCorr_Table,['dataSet' num2str(dataSet) '_networkMeasures_comparisons.xlsx']);

networkMeasures_Polychoric = networkMeasures;
networkMeasureMeans_Polychoric = networkMeasureMeans;
networkMeasures_Pearson = networkMeasures;
networkMeasureMeans_Pearson = networkMeasureMeans;
networkMeasureWeightedMeans_Pearson = networkMeasureWeightedMeans;
networkMeasures_MI = networkMeasures;
networkMeasureMeans_MI = networkMeasureMeans;
networkMeasureWeightedMeans_MI = networkMeasureWeightedMeans;

for i_var = 1:size(networkMeasures,2)
    netCorr(1,i_var) = corr(networkMeasures_Pearson(:,i_var),networkMeasures_Polychoric(:,i_var));
    netCorr(2,i_var) = corr(networkMeasures_Pearson(:,i_var),networkMeasures_MI(:,i_var));
    netCorr(3,i_var) = corr(networkMeasures_Polychoric(:,i_var),networkMeasures_MI(:,i_var));
    netCorr(4,i_var) = corr(networkMeasureMeans_Pearson(:,i_var),networkMeasureMeans_Polychoric(:,i_var));
    netCorr(5,i_var) = corr(networkMeasureMeans_Pearson(:,i_var),networkMeasureMeans_MI(:,i_var));
    netCorr(6,i_var) = corr(networkMeasureMeans_Polychoric(:,i_var),networkMeasureMeans_MI(:,i_var));
    netCorr(7,i_var) = corr(networkMeasureWeightedMeans_Pearson(:,i_var),networkMeasureWeightedMeans_MI(:,i_var));
end
rowNames2 = {'staticR v staticP','staticR v staticMI','staticP v staticMI','aggR v aggP','aggR v aggMI','aggP v aggMI','weightR v weightMI'}';
netCorr_Table = array2table(netCorr,'VariableNames',variableNames);
netCorr_Table = horzcat(rowNames2,netCorr_Table);
netCorr_Table.Properties.VariableNames{1} = 'comparison';
writetable(netCorr_Table,['dataSet' num2str(dataSet) '_networkMeasures_comparisons_across types.xlsx']);