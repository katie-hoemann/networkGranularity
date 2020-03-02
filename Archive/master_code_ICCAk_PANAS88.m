clear;
clc;

filename = '88PANAS_RAW.dat';
rawdata = importdata(filename);
subjectIDList = unique(rawdata.data(:,1));

load('wordslist.mat', 'words','M')
ICC_mat = zeros(length(subjectIDList),4);

%{
[A, I1] = sort(words);
[B, I2] = sort(M(:,1));
[~, I3] = sort(I1);

strcmp(A, B);
M = M(I2,:);
M = M(I3,:);
%}

labels = cellfun(@conv_label,M(:,2:end));

%%
for i_subject = 1:length(subjectIDList)
    data = [];
    subjectID = subjectIDList(i_subject);
    index = find(rawdata.data(:,1)==subjectID);
    data = rawdata.data(index,2:end);
    score = data(:,3:end);
    %% Removing missing data
    missing_data = isnan(score);
    missing_data2 = any(missing_data,2);
    score = score(~missing_data2,:);
    
    %% Computation of ICC
    ICC_mat(i_subject,1) = ICC(score(:,(labels(:,1)==0)),'A-k');
    ICC_mat(i_subject,2) = ICC(score(:,(labels(:,1)==1)),'A-k');
    ICC_mat(i_subject,3) = ICC(score(:,(labels(:,2)==0)),'A-k');
    ICC_mat(i_subject,4) = ICC(score(:,(labels(:,2)==1)),'A-k');
end

ReportTable = array2table(ICC_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
    'Pos_V','Low_A','High_A'});
writetable(ReportTable,'ICC_Ak.dat','WriteRowNames',true);
writetable(ReportTable,'ICC_Ak.txt','WriteRowNames',true);

function y = conv_label(x)
switch x
    case {'Negative','Low'}
        y = 0;
    case {'Positive','High'}
        y = 1;
    otherwise
        y = 4;
end
end

