clear all;
clc;

%% specify dataset and parameters
dataSet = 39; % must be set to 18
SDs = 0; % set to 1 to use standard deviations
fontsize = 15; % set fontsize for axis labels

%% import data
if dataSet == 18
    data = importdata('18term_Pearson_windowed.xlsx');
    ZmGAD7 = data.data(:,148); % standardized GAD7 scores
    ZmPHQ8 = data.data(:,149); % standardized PHQ8 scores 
    moodSymptoms = (ZmGAD7+ZmPHQ8)/2;
    if SDs == 1
        granFactor = data.data(:,156);
    else
        granFactor = data.data(:,155);
    end
elseif dataSet == 39
    data = importdata('39term_Pearson_windowed.xlsx');
    ZBAI = data.data(:,134); % standardized BAI scores
    moodSymptoms = ZBAI;
    if SDs == 1
        granFactor = data.data(:,144);
    else
        granFactor = data.data(:,140); % use factor #2 but label #3
    end
end

%% create scatter plots of significant relationships
if dataSet == 18
    if SDs == 1
        figure;
        scatter1 = scatter(granFactor,moodSymptoms,[],rgb('MediumPurple'),'filled');
        xlabel('granularity factor #1 (SD estimates)');
        xlim([-3 3]);
        ylim([-3 3]);
        set(gca,'fontsize',fontsize)
        ylabel('self-reported mood symptoms');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_scatter_plot_window_SDs'];
        saveas(scatter1,filename,'tiff');
    else
        figure;
        scatter1 = scatter(granFactor,moodSymptoms,[],rgb('MediumPurple'),'filled');
        xlabel('granularity factor #2 (mean estimates)');
        xlim([-3 3]);
        ylim([-3 3]);
        set(gca,'fontsize',fontsize)
        ylabel('self-reported mood symptoms');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_scatter_plot_window_Ms'];
        saveas(scatter1,filename,'tiff');
    end
elseif dataSet == 39
    if SDs == 1
        figure;
        scatter1 = scatter(granFactor,moodSymptoms,[],rgb('SeaGreen'),'filled');
        xlabel('granularity factor (SD estimates)');
        xlim([-3 3]);
        ylim([-3 3]);
        set(gca,'fontsize',fontsize)
        ylabel('self-reported mood symptoms');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_scatter_plot_window_SDs'];
        saveas(scatter1,filename,'tiff');
    else
        figure;
        scatter1 = scatter(granFactor,moodSymptoms,[],rgb('SeaGreen'),'filled');
        xlabel('granularity factor #3 (mean estimates)');
        xlim([-3 3]);
        ylim([-3 3]);
        set(gca,'fontsize',fontsize)
        ylabel('self-reported mood symptoms');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_scatter_plot_window_Ms'];
        saveas(scatter1,filename,'tiff');
    end
end

close all;
