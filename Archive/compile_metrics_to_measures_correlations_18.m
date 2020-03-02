clear all;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
print = 0; % set to 1 to generate spreadsheets with results

%% load data
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

rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASItotal = rawData.data(:,7)-21; % subtract 21 to bring in line with original rating scale
ERStotal = rawData.data(:,11); 
GADtotal = rawData.data(:,12); % 12 for session 1; 34 for session 2
PSStotal = rawData.data(:,21).*3.5; %46 for session 2 % multiply by 3.5 to bring in line with norms from 14-item measure
PHQ15total = rawData.data(:,22); %47 for session 2
PHQ8total = rawData.data(:,23); %48 for session 2
neuroticism = rawData.data(:,49);
TAStotal = rawData.data(:,17); %39 for session 2
RDEES = rawData.data(:,20)./14; %45 for session 2 % divide by 14 to get mean per item

%% network metrics of interest
x1 = numCommunities(:,2); % number of communities at no threshold
x2 = modularity_VA(:,2); % modularity for valence-assigned communities at no threshold
x3 = clusterCoefGlobal(:,2); % global clustering at no threshold
x4 = networkDensity(:,5); % weighted density at no threshold
x5 = diversity(:,2); % positive diversity at no threshold
x6 = diversity(:,5); % negative diversity at no threshold

metrics = [x1 x2 x3 x4 x5 x6];
metricNames = {'numCom' 'modVA' 'cluster' 'density' 'posDiv' 'negDiv'};
[~,numMetrics] = size(metrics);

%% psychological measures 
y1 = ERStotal(:);
y2 = PHQ8total(:);
y3 = PSStotal(:).*3.5; % multiply by 3.5 to bring in line with 14-item measure
y4 = GADtotal(:);
y5 = neuroticism(:);
y6 = PHQ15total(:);
y7 = ASItotal(:);
y8 = TAStotal(:);
y9 = RDEES(:);

measures = [y1 y2 y3 y4 y5 y6 y7 y8 y9];
measureNames = {'ERS' 'PHQ8' 'PSS' 'GAD' 'neuro' 'PHQ15' 'ASI' 'TAS' 'RDEES'};
[~,numMeasures] = size(measures);

%% other variables to test
z1 = rtoZ(:,3); % zICC total
z2 = affect(:,2); % mean negative affect
z3 = affect(:,4); % std negative affect
z4 = affect(:,1); % mean positive affect
z5 = affect(:,3); % std positive affect

allPredictors = [metrics z1 z2 z3 z4 z5];
predictorNames = {'numCom' 'modVA' 'cluster' 'density' 'posDiv' 'negDiv' 'zICC' 'mNeg' 'sdNeg' 'mPos' 'sdPos'};

%% compile all raw data and export
allData = [allPredictors measures];
allDataTable = array2table(allData);
allDataTable.Properties.VariableNames = [predictorNames measureNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_all_raw_data.xlsx'];
    writetable(allDataTable,filename);
end

%% create correlation plots
for i_measure = 1:numMeasures
    y = measures(:,i_measure);   
    toPlot = [y allPredictors];
    toPlot = array2table(toPlot);
    toPlot.Properties.VariableNames = [measureNames{i_measure} predictorNames];
    figure;
    corrplot(toPlot)
end

%% compute zero-order correlation matrices
for i_measure = 1:numMeasures
    toRun = [measures(:,i_measure) metrics];
    [corrs corrPs] = corrcoef(toRun);
    corrMatrices(:,:,i_measure) = corrs;
    corrsTable = array2table(corrs);
    corrsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_zero_order_correlations.xlsx'];
        writetable(corrsTable,filename);
    end
    corrPsMatrices(:,:,i_measure) = corrPs;
    corrPsTable = array2table(corrPs);
    corrPsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_zero_order_correlations_pvalues.xlsx'];
        writetable(corrPsTable,filename);
    end
end

%% compute partial correlation matrices
for i_measure = 1:numMeasures
    toRun = [measures(:,i_measure) metrics];
    [partials partialPs] = partialcorr(toRun);
    partialsMatrices(:,:,i_measure) = partials;
    partialsTable = array2table(partials);
    partialsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_partial_correlations.xlsx'];
        writetable(partialsTable,filename);
    end
    partialPsMatrices(:,:,i_measure) = partialPs;
    partialPsTable = array2table(partialPs);
    partialPsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_partial_correlations_pvalues.xlsx'];
        writetable(partialPsTable,filename);
    end
end

zICCtoRun = [z1 metrics];
[zICCpartials zICCpartialPs] = partialcorr(zICCtoRun);
zICCpartialsTable = array2table(zICCpartials);
zICCpartialsTable.Properties.VariableNames = ['zICC' metricNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_zICC_partial_correlations.xlsx'];
    writetable(zICCpartialsTable,filename);
end
zICCpartialPsTable = array2table(zICCpartialPs);
zICCpartialPsTable.Properties.VariableNames = ['zICC' metricNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_zICC_partial_correlations_pvalues.xlsx'];
    writetable(zICCpartialPsTable,filename);
end

%% run regressions (to get residualized variables)
% a partial correlation is computed between two residuals; a semipartial
% correlation is computed between one residual and one raw variable
% http://faculty.cas.usf.edu/mbrannick/regression/Partial.html
% note this means the residuals need to take into account the same IVs

% for i_measure = 1:numMeasures
%     y = measures(:,i_measure);
%     X = [ones(length(y),1) metrics];
%     [bY,~,residY,~,statsY] = regress(y,X);
%     regressCoefY(:,i_measure) = bY;
%     regressResidY(:,i_measure) = residY;
%     regressStatsY(:,i_measure) = statsY';
% end

for i_metric = 1:numMetrics
    x = metrics(:,i_metric);
    x_copy = metrics;
    x_copy(:,i_metric) = [];
    X = [ones(length(x),1) x_copy];
    [bX,~,residX,~,statsX] = regress(x,X);
    regressCoefX(:,i_metric) = bX;
    regressResidX(:,i_metric) = residX;
    regressStatsX(:,i_metric) = statsX';
end

%% compute semi-partial correlation matrices (using residualized metrics)
for i_measure = 1:numMeasures
    toRun = [measures(:,i_measure) regressResidX];
    [semipartials semiPs] = corrcoef(toRun);
    semipartialMatrices(:,:,i_measure) = semipartials;
    semipartialsTable = array2table(semipartials);
    semipartialsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_semi-partial_correlations.xlsx'];
        writetable(semipartialsTable,filename);
    end
    semiPsMatrices(:,:,i_measure) = semiPs;
    semiPsTable = array2table(semiPs);
    semiPsTable.Properties.VariableNames = [measureNames{i_measure} metricNames];
    if print == 1
        filename = ['dataSet' num2str(dataSet) '_' char(measureNames{i_measure}) '_semi-partial_correlations_pvalues.xlsx'];
        writetable(semiPsTable,filename);
    end
end

zICCtoRun = [z1 regressResidX];
[zICCsemipartials zICCsemiPs] = corrcoef(zICCtoRun);
zICCsemipartialsTable = array2table(zICCsemipartials);
zICCsemipartialsTable.Properties.VariableNames = ['zICC' metricNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_zICC_semi-partial_correlations.xlsx'];
    writetable(zICCsemipartialsTable,filename);
end
zICCsemiPsTable = array2table(zICCsemiPs);
zICCsemiPsTable.Properties.VariableNames = ['zICC' metricNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_zICC_semi-partial_correlations_pvalues.xlsx'];
    writetable(zICCsemiPsTable,filename);
end
