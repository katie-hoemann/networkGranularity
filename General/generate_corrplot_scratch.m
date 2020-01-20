load('dataSet18_networkMeasures.mat')
toPlot = networkMeasureWeightedSDs;
[2 5 6 7 8 9 12 13];
toPlot(:,ans) = [];
toPlot(:,5) = [];
%corrplot(toPlot,'varNames',{'cluster','density','numCom','posDiv','negDiv','comRadius','modVA'})
corrplot(toPlot,'varNames',{'cluster','density','numCom','posDiv','comRadius','modVA'})