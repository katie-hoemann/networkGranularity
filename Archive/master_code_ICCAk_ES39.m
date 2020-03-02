clear;
clc;

filename = 'ES39data.mat';
rawdata = load(filename); %edit command b/c loading .mat instead of .dat
subjectIDList = unique(rawdata.ESDATA(:,1)); %edit substructure name

%load('wordslist.mat', 'words','M') %remove; new word list file loaded next
words39=readtable('words39.csv'); %load word list, valence and arousal values

for i=1:height(words39)
if words39.Valence(i)>5
    ValCat(i)={'Positive'};
    Positive(i)=1;
else
    ValCat(i)={'Negative'};
    Positive(i)=0;
end
end %define valence categories

for i=1:height(words39)
if words39.Arousal(i)>4.6
    AroCat(i)={'High'};
    High(i)=1;
else
    AroCat(i)={'Low'};
    High(i)=0;
end
end %define arousal categories

words39=[words39 ValCat' AroCat']; %append table
words39.Properties.VariableNames(5:end)={'ValCat' 'AroCat'}; %label new vars
labels=[Positive' High']; %create matrix for logical indexing in ICC commands

rawICC_mat = zeros(length(subjectIDList),8);
gran_mat = zeros(length(subjectIDList),8);
rtoz_mat = zeros(length(subjectIDList),8);
zinv_mat = zeros(length(subjectIDList),8);

%{
[A, I1] = sort(words);
[B, I2] = sort(M(:,1));
[~, I3] = sort(I1);

strcmp(A, B);
M = M(I2,:);
M = M(I3,:);
%} 

%labels = cellfun(@conv_label,words39(:,5:end));  %remove; replaced above

%%
for i_subject = 1:length(subjectIDList)
    data = [];
    subjectID = subjectIDList(i_subject);
    index = find(rawdata.ESDATA(:,1)==subjectID);
    data = rawdata.ESDATA(index,2:end);
    %score = data(:,:); %remove; unnecessary given updated source data file
    %% Removing missing data
    missing_data = isnan(data); %edit 'score' to 'data'
    missing_data2 = any(missing_data,2); %edit 'score' to 'data'
    data = data(~missing_data2,:); %edit 'score' to 'data'
    %% Rescale to from 1-5 to 0-4 %should add an if clause so only do this if scale doesn't start at 0
    data = data-1;
    %% Computation of ICCs - by hemispheres and quadrants
    rawICC_mat(i_subject,1) = ICC(data(:,(labels(:,1)==0)),'A-k'); %edit 'score' to 'data'
    rawICC_mat(i_subject,2) = ICC(data(:,(labels(:,1)==1)),'A-k'); %edit 'score' to 'data'
    rawICC_mat(i_subject,3) = ICC(data(:,(labels(:,2)==0)),'A-k'); %edit 'score' to 'data'
    rawICC_mat(i_subject,4) = ICC(data(:,(labels(:,2)==1)),'A-k'); %edit 'score' to 'data'
    rawICC_mat(i_subject,5) = ICC(data(:,(labels(:,1)==0 & labels(:,2)==0)),'A-k');
    rawICC_mat(i_subject,6) = ICC(data(:,(labels(:,1)==1 & labels(:,2)==1)),'A-k');
    rawICC_mat(i_subject,7) = ICC(data(:,(labels(:,1)==0 & labels(:,2)==1)),'A-k');
    rawICC_mat(i_subject,8) = ICC(data(:,(labels(:,1)==1 & labels(:,2)==0)),'A-k');
    %% Calculation of granularity from ICCs
    gran_mat(i_subject,1) = 1-rawICC_mat(i_subject,1);
    gran_mat(i_subject,2) = 1-rawICC_mat(i_subject,2);
    gran_mat(i_subject,3) = 1-rawICC_mat(i_subject,3);
    gran_mat(i_subject,4) = 1-rawICC_mat(i_subject,4);
    gran_mat(i_subject,5) = 1-rawICC_mat(i_subject,5);
    gran_mat(i_subject,6) = 1-rawICC_mat(i_subject,6);
    gran_mat(i_subject,7) = 1-rawICC_mat(i_subject,7);
    gran_mat(i_subject,8) = 1-rawICC_mat(i_subject,8);
    %% Fisher transform ICCs
    rtoz_mat(i_subject,1) = 0.5*log((1+rawICC_mat(i_subject,1))/(1-rawICC_mat(i_subject,1)));
    rtoz_mat(i_subject,2) = 0.5*log((1+rawICC_mat(i_subject,2))/(1-rawICC_mat(i_subject,2)));
    rtoz_mat(i_subject,3) = 0.5*log((1+rawICC_mat(i_subject,3))/(1-rawICC_mat(i_subject,3)));
    rtoz_mat(i_subject,4) = 0.5*log((1+rawICC_mat(i_subject,4))/(1-rawICC_mat(i_subject,4)));
    rtoz_mat(i_subject,5) = 0.5*log((1+rawICC_mat(i_subject,5))/(1-rawICC_mat(i_subject,5)));
    rtoz_mat(i_subject,6) = 0.5*log((1+rawICC_mat(i_subject,6))/(1-rawICC_mat(i_subject,6)));
    rtoz_mat(i_subject,7) = 0.5*log((1+rawICC_mat(i_subject,7))/(1-rawICC_mat(i_subject,7)));
    rtoz_mat(i_subject,8) = 0.5*log((1+rawICC_mat(i_subject,8))/(1-rawICC_mat(i_subject,8)));
    %% Invert Fisher-transformed ICCs for intuitive directionality
    zinv_mat(i_subject,1) = rtoz_mat(i_subject,1)*-1;
    zinv_mat(i_subject,2) = rtoz_mat(i_subject,2)*-1;
    zinv_mat(i_subject,3) = rtoz_mat(i_subject,3)*-1;
    zinv_mat(i_subject,4) = rtoz_mat(i_subject,4)*-1;
    zinv_mat(i_subject,5) = rtoz_mat(i_subject,5)*-1;
    zinv_mat(i_subject,6) = rtoz_mat(i_subject,6)*-1;
    zinv_mat(i_subject,7) = rtoz_mat(i_subject,7)*-1;
    zinv_mat(i_subject,8) = rtoz_mat(i_subject,8)*-1;
end

ReportTableICC = array2table(rawICC_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
    'Pos_V','Low_A','High_A','NegLow','PosHigh','NegHigh','PosLow'});
writetable(ReportTableICC,'ICC_Ak.dat','WriteRowNames',true);
writetable(ReportTableICC,'ICC_Ak.txt','WriteRowNames',true);

ReportTableGran = array2table(gran_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
    'Pos_V','Low_A','High_A','NegLow','PosHigh','NegHigh','PosLow'});
writetable(ReportTableGran,'Gran_Ak.dat','WriteRowNames',true);
writetable(ReportTableGran,'Gran_Ak.txt','WriteRowNames',true);

ReportTableRtoZ = array2table(rtoz_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
    'Pos_V','Low_A','High_A','NegLow','PosHigh','NegHigh','PosLow'});
writetable(ReportTableRtoZ,'RtoZ_Ak.dat','WriteRowNames',true);
writetable(ReportTableRtoZ,'RtoZ_Ak.txt','WriteRowNames',true);

ReportTableZInv = array2table(zinv_mat,'RowNames',cellstr(num2str(subjectIDList)),'VariableNames',{'Neg_V',...
    'Pos_V','Low_A','High_A','NegLow','PosHigh','NegHigh','PosLow'});
writetable(ReportTableZInv,'ZInv_Ak.dat','WriteRowNames',true);
writetable(ReportTableZInv,'ZInv_Ak.txt','WriteRowNames',true);

%{
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
%}

