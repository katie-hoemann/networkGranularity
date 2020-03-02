clear all;
clc;

%% specify dataset
dataSet = 88; 
print = 1; % set to 1 to generate spreadsheets with results

%% load data
filenames{1} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_instability.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_NookMeasures.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end

if dataSet == 39
    emodiversity(29,:) = [];
    instability(29,:) = [];
    NookMeasures(29,:) = [];
end

%% compile all raw data and export
allData = [emodiversity(:) instability(:) NookMeasures(:,1) NookMeasures(:,2) NookMeasures(:,3)];
allDataTable = array2table(allData);
allDataTable.Properties.VariableNames = {'emodiversity' 'instability' 'singleEmotion' 'intensity' 'extremeScaleUse'};
if print == 1
    filename = ['dataSet' num2str(dataSet) '_additional_derived_measures.xlsx'];
    writetable(allDataTable,filename);
end