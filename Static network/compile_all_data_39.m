clear all;
clc;

%% specify dataset
dataSet = 39; % must be set to 18
print = 1; % set to 1 to generate spreadsheets with results

%% load data
filenames{1} = ['dataSet' num2str(dataSet) '_communities.mat'];
filenames{2} = ['dataSet' num2str(dataSet) '_modularity.mat'];
filenames{3} = ['dataSet' num2str(dataSet) '_modularity_VA.mat'];
filenames{4} = ['dataSet' num2str(dataSet) '_clustering.mat'];
filenames{5} = ['dataSet' num2str(dataSet) '_clustering_SD.mat'];
filenames{6} = ['dataSet' num2str(dataSet) '_participation.mat'];
filenames{7} = ['dataSet' num2str(dataSet) '_participation_SD.mat'];
filenames{8} = ['dataSet' num2str(dataSet) '_gateway.mat'];
filenames{9} = ['dataSet' num2str(dataSet) '_gateway_SD.mat'];
filenames{10} = ['dataSet' num2str(dataSet) '_diversity.mat'];
filenames{11} = ['dataSet' num2str(dataSet) '_diversity_SD.mat'];
filenames{12} = ['dataSet' num2str(dataSet) '_density.mat'];
filenames{13} = ['dataSet' num2str(dataSet) '_comm_radius.mat'];
filenames{14} = ['dataSet' num2str(dataSet) '_zICC.mat'];
filenames{15} = ['dataSet' num2str(dataSet) '_affect.mat'];
filenames{16} = ['dataSet' num2str(dataSet) '_emodiversity.mat'];
filenames{17} = ['dataSet' num2str(dataSet) '_instability.mat'];
filenames{18} = ['dataSet' num2str(dataSet) '_NookMeasures.mat'];
for i_file = 1:numel(filenames)
    load(filenames{i_file});
end

rawData = importdata('39ESdata_QuestionnaireData.xlsx');
BAItotal = rawData.data(:,9); 
BDItotal = rawData.data(:,12);
SWLStotal = rawData.data(:,13);
TAStotal = rawData.data(:,2);
TASddf = rawData.data(:,3);
TASdif = rawData.data(:,4);
TASeot = rawData.data(:,5);
RDEES = rawData.data(:,6);
diff = rawData.data(:,7);
range = rawData.data(:,8);

%% remove subject 31 (row 29)
numCommunities(29,:) = [];
modularity(29,:) = [];
modularity_VA(29,:) = [];
clusterCoefGlobal(29,:) = [];
clusterCoefGlobalSD(29,:) = [];
participation(29,:) = [];
participationSD(29,:) = [];
gateway(29,:) = [];
gatewaySD(29,:) = [];
diversity(29,:) = [];
diversitySD(29,:) = [];
networkDensity(29,:) = [];
commRadius(29,:) = [];
rtoZ(29,:) = [];
affect(29,:) = [];
emodiversity(29,:) = [];
instability(29,:) = [];
NookMeasures(29,:) = [];
BAItotal(29,:) = [];
BDItotal(29,:) = [];
SWLStotal(29,:) = [];
TAStotal(29,:) = [];
TASddf(29,:) = [];
TASdif(29,:) = [];
TASeot(29,:) = [];
RDEES(29,:) = [];
diff(29,:) = [];
range(29,:) = [];

%% network metrics of interest
x1 = numCommunities(:,2); % number of communities at no threshold
x2 = modularity_VA(:,2); % modularity for valence-assigned communities at no threshold
x3 = clusterCoefGlobal(:,2); % global clustering at no threshold
x4 = clusterCoefGlobalSD(:,2); % SD clustering at no threshold
x5 = networkDensity(:,5); % weighted density at no threshold
x6 = diversity(:,2); % positive diversity at no threshold
x7 = diversitySD(:,2); % SD positive diversity at no threshold
x8 = diversity(:,5); % negative diversity at no threshold
x9 = diversitySD(:,5); % SD negative diversity at no threshold
x10 = commRadius(:,2); % average community radius

networkMetrics = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10];
networkMetricNames = {'numCom' 'modVA' 'cluster' 'clusterSD' 'density' 'posDiv' 'posDivSD' 'negDiv' 'negDivSD' 'comRadius'};
[~,numMetrics] = size(networkMetrics);

%% psychological measures 
y1 = BAItotal(:);
y2 = BDItotal(:);
y3 = SWLStotal(:);
y4 = TAStotal(:);
y5 = TASddf(:);
y6 = TASdif(:);
y7 = TASeot(:);
y8 = RDEES(:);
y9 = diff(:);
y10 = range(:);

psychMeasures = [y1 y2 y3 y4 y5 y6 y7 y8 y9 y10];
psychMeasureNames = {'BAI' 'BDI' 'SWLS' 'TAS' 'DDF' 'DIF' 'EOT' 'RDEES' 'diff' 'range'};
[~,numMeasures] = size(psychMeasures);

%% other derived measures
z1 = rtoZ(:,3); % zICC total
z2 = affect(:,2); % mean negative affect
z3 = affect(:,4); % std negative affect
z4 = affect(:,1); % mean positive affect
z5 = affect(:,3); % std positive affect
z6 = emodiversity(:); % emodiversity (Quoidbach formula)
z7 = instability(:); % instability (from Dejonkheere et al)
z8 = NookMeasures(:,1); % single emotion rating (from Nook et al)
z9 = NookMeasures(:,2); % average intensity (from Nook et al)
z10 = NookMeasures(:,3); % extreme scale use (from Nook et al)

derivedMeasures = [z1 z2 z3 z4 z5 z6 z7 z8 z9 z10];
derivedMeasureNames = {'zICC' 'mNeg' 'sdNeg' 'mPos' 'sdPos' 'emodiversity' 'instability' 'singleEmotion' 'intensity' 'extremeScaleUse'};
[~,numDerived] = size(derivedMeasures);

%% compile all raw data and export
allData = [networkMetrics derivedMeasures psychMeasures];
allDataTable = array2table(allData);
allDataTable.Properties.VariableNames = [networkMetricNames derivedMeasureNames psychMeasureNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_all_raw_data.xlsx'];
    writetable(allDataTable,filename);
end

%% run Kolmogorov-Smirnov test to assess normality; normalize distribution as needed
numVariables = numMetrics + numMeasures + numDerived;
for i_variable = 1:numVariables
    [h(i_variable),p(i_variable),ks(i_variable)] = lillietest(allData(:,i_variable));
    if h(i_variable) == 1
        fracRank = FractionalRankings(allData(:,i_variable));
        fracRank = fracRank/(max(fracRank));
        fracRank(fracRank == 1) = .999;
        fracRankNorm(:,i_variable) = norminv(fracRank,mean(allData(:,i_variable)),std(allData(:,i_variable)));
    elseif h(i_variable) == 0
        fracRankNorm(:,i_variable) = NaN(length(allData),1);
    end
end

%% compile all normalized data and export
fracRankNormTable = array2table(fracRankNorm);
fracRankNormTable.Properties.VariableNames = [networkMetricNames derivedMeasureNames psychMeasureNames];
if print == 1
    filename = ['dataSet' num2str(dataSet) '_all_normalized_data.xlsx'];
    writetable(fracRankNormTable,filename);
end