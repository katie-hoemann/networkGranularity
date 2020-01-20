function [R2change, Fchange, pChange] = hmlr(y, xStep1, xStep2)

% https://www.ibm.com/support/knowledgecenter/en/SSLVMB_22.0.0/com.ibm.spss.statistics.algorithms/alg_regression_fchange.htm
% https://www.researchgate.net/post/What_is_a_significant_f_change_value_in_a_hierarchical_multiple_regression

X_padded = [ones(length(y),1) xStep1 xStep2];
X_step1 = X_padded(:,1:2);
X_step2 = X_padded;
idx = find(isnan(X_padded));
[row,~] = ind2sub(size(X_padded),idx);
X_step1(row,:) = [];
X_step2(row,:) = [];
y(row) = [];
[~,~,~,~,stats_step1] = regress(y,X_step1);
[~,~,~,~,stats_step2] = regress(y,X_step2);

R2_step1 = stats_step1(1);
R2_step2 = stats_step2(1);
[~,IVs_step1] = size(X_step1);
[~,IVs_step2] = size(X_step2);
IVs_step1 = IVs_step1-1;
IVs_step2 = IVs_step2-1;
df1 = IVs_step2-IVs_step1; % number of IVs added
df2 = length(X_padded)-IVs_step2-1; % N-k-1
R2change = R2_step2-R2_step1;
Fchange = [(R2_step2-R2_step1)/df1]/[(1-R2_step2)/df2];
pChange = 1-fcdf(F_change,df1,df2);