x1 = rtoZ(:,3); % zICCtotal
x2 = affect(:,2); % mean negative affect
x3 = affect(:,4); % std negative affect
x4 = affect(:,1); % mean positive affect
x5 = affect(:,3); % std positive affect
x6 = diversity(:,5); % negative diversity at no threshold
x7 = modularity_VA(:,2); % modularity for valence-assigned communities at no threshold
x8 = networkDensity(:,5); % weighted network density at no threshold
x9 = participation(:,2); % positive participation at no threshold
x10 = modularity(:,2); % modularity for detected communities at no threshold
x11 = numCommunities(:,2); % number of communities at no threshold
x12 = diversity(:,2); % positive diversity at no threshold
y1 = ERStotal(:);
y2 = PHQ8total(:);
y3 = PSStotal(:).*3.5; % multiply by 3.5 to bring in line with 14-item measure
y4 = GADtotal(:);
y5 = neuroticism(:);
y6 = PHQ15total(:);

%% hierarchical multiple regression
% https://www.ibm.com/support/knowledgecenter/en/SSLVMB_22.0.0/com.ibm.spss.statistics.algorithms/alg_regression_fchange.htm
% https://www.researchgate.net/post/What_is_a_significant_f_change_value_in_a_hierarchical_multiple_regression

X_padded = [ones(length(y3),1) x1 x10 x9];
X_step1 = X_padded(:,1:2);
X_step2 = X_padded;
idx = find(isnan(X_padded));
[row,~] = ind2sub(size(X_padded),idx);
X_step1(row,:) = [];
X_step2(row,:) = [];
y3(row) = [];
[b_step1,~,r_step1,~,stats_step1] = regress(y3,X_step1);
[b_step2,~,r_step2,~,stats_step2] = regress(y3,X_step2);

R2_step1 = stats_step1(1);
R2_step2 = stats_step2(1);
[~,IVs_step1] = size(X_step1);
[~,IVs_step2] = size(X_step2);
IVs_step1 = IVs_step1-1;
IVs_step2 = IVs_step2-1;
df1 = IVs_step2-IVs_step1; % number of IVs added
df2 = length(X_padded)-IVs_step2-1; % N-k-1
R2_change = R2_step2-R2_step1
F_change = [(R2_step2-R2_step1)/df1]/[(1-R2_step2)/df2]
p_change = 1-fcdf(F_change,df1,df2)

%% stepwise regression
% X = [x1 x6 x7];
% data1 = [x1 x6 y1];
% data1 = array2table(data1);
% data1.Properties.VariableNames = {'zICC','negDiversity','ERStotal'};
% data2 = [x1 x7 y2];
% data2 = array2table(data2);
% data2.Properties.VariableNames = {'zICC','modularity_VA','PHQ8total'};
% data3 = [x1 x7 y3];
% data3 = array2table(data3);
% data3.Properties.VariableNames = {'zICC','modularity_VA','PSStotal'};
% 
% model1 = stepwiselm(data3,'constant','ResponseVar','PSStotal')
% model2 = stepwiselm(data3,'interactions','ResponseVar','PSStotal')
% 
% [b,~,pval,inmodel,stats,~,history] = stepwisefit(X,y3)