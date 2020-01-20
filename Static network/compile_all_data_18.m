clear all;
clc;

%% specify dataset
dataSet = 18; % must be set to 18
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

rawData = importdata('18ARIVariations_QuestionnaireData.xlsx');
ASItotal = rawData.data(:,7)-21; % subtract 21 to bring in line with original rating scale
ERStotal = rawData.data(:,11); 
GADtotal = rawData.data(:,12); 
PSStotal = rawData.data(:,21).*3.5; %21 for session 1; 46 for session 2 % multiply by 3.5 to bring in line with norms from 14-item measure
PHQ15total = rawData.data(:,22); %22 for session 1; 47 for session 2
PHQ8total = rawData.data(:,23); %23 for session 1; 48 for session 2
neuroticism = rawData.data(:,49);
TAStotal = rawData.data(:,39); %17 for session 1; 39 for session 2
TASeot = rawData.data(:,36); %14 for session 1; 36 for session 2
TASddf = rawData.data(:,37); %15 for session 1; 37 for session 2
TASdif = rawData.data(:,38); %16 for session 1; 38 for session 2
RDEES = rawData.data(:,45)./14; %20 for session 1; 45 for session 2 % divide by 14 to get mean per item
diff = rawData.data(:,43)./14; %18 for session 1; 43 for session 2 % divide by 14 to get mean per item
range = rawData.data(:,44)./14; %19 for session 1; 44 for session 2 % divide by 14 to get mean per item

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
y1 = ERStotal(:);
y2 = PHQ8total(:);
y3 = PSStotal(:).*3.5; % multiply by 3.5 to bring in line with 14-item measure
y4 = GADtotal(:);
y5 = neuroticism(:);
y6 = PHQ15total(:);
y7 = ASItotal(:);
y8 = TAStotal(:);
y9 = TASeot(:);
y10 = TASddf(:);
y11 = TASdif(:);
y12 = RDEES(:);
y13 = diff(:);
y14 = range(:);

psychMeasures = [y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12 y13 y14];
psychMeasureNames = {'ERS' 'PHQ8' 'PSS' 'GAD' 'neuro' 'PHQ15' 'ASI' 'TAS' 'EOT' 'DDF' 'DIF' 'RDEES' 'diff' 'range'};
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