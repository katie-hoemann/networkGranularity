clear all;
clc;

%% specify dataset
dataSet = 39; % must be set to 39
print = 1; % set to 1 to generate spreadsheets with results

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

%% load criterion validity measures
rawData = importdata('39ESdata_QuestionnaireData.xlsx');
BAItotal = rawData.data(:,9); 
BDItotal = rawData.data(:,12);
SWLStotal = rawData.data(:,13);
TAStotal = rawData.data(:,2);
RDEES = rawData.data(:,6);

%% remove subject 31 (row 29) for all calculations
rtoZ(29,:) = [];
numCommunities(29,:) = [];
modularity(29,:) = [];
modularity_VA(29,:) = [];
clusterCoefGlobal(29,:) = [];
participation(29,:) = [];
gateway(29,:) = [];
diversity(29,:) = [];
networkDensity(29,:) = [];
commRadius(29,:) = [];
affect(29,:) = [];
BAItotal(29,:) = [];
BDItotal(29,:) = [];
SWLStotal(29,:) = [];
TAStotal(29,:) = [];
RDEES(29,:) = [];

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
y1 = BAItotal(:);
y2 = BDItotal(:);
y3 = SWLStotal(:);
y4 = TAStotal(:);
y5 = RDEES(:);

measures = [y1 y2 y3 y4 y5];
measureNames = {'BAI' 'BDI' 'SWLS' 'TAS' 'RDEES'};
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
