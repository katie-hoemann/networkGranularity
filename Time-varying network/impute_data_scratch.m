test = gran_M_array(2:end,2:end);
deletedPPs = zeros(size(test,1),1);
for i_row = 1:size(test,1)
    testRow = test(i_row,:);
    testRow = testRow(~isnan(testRow));
    if numel(testRow) < 8
        deletedPPs(i_row) = 1;
    end
end
deletedPPs = find(deletedPPs==1);
test(deletedPPs,:) = [];

fillMean = fillmissing(test,'movmean',3,2);
%fillLinear = fillmissing(test,'linear',2,'EndValues','nearest');

moveMean = movemean(test,3,2,'omitnan','Endpoints','discard');
%this is more for smoothing the data rather than calculating windows

%weighting for days: day1*1/6, day2*4/6, day3*1/6