clear all;
clc;

%% specify dataset and parameters
dataSet = 39; % must be set to 18
windowed = 0; % set to 1 to use data from sliding windows
SDs = 0; % set to 1 to use standard deviations from windowed data
polychoric = 1; % set to 1 to use polychoric networks (static only)
fontsize = 15; % set fontsize for axis labels

%% import data
if dataSet == 18
    if windowed == 1
        data = importdata('18term_Pearson_windowed.xlsx');
        mGAD7 = data.data(:,62);
        mPHQ8 = data.data(:,58);
        if SDs == 1
            granFactor = data.data(:,156);
        else
            granFactor = data.data(:,155);
        end
    elseif polychoric == 1
        data = importdata('18term_polychoric_static.xlsx');
        mGAD7 = data.data(:,28);
        mPHQ8 = data.data(:,24);
        granFactor = data.data(:,87);
    else
        data = importdata('18term_Pearson_static.xlsx');
        mGAD7 = data.data(:,28);
        mPHQ8 = data.data(:,24);
        granFactor = data.data(:,98);
    end
elseif dataSet == 39
    if windowed == 1
        data = importdata('39term_Pearson_windowed.xlsx');
        BAI = data.data(:,55);
        BDI = data.data(:,56);
        if SDs == 1
            granFactor = data.data(:,144);
        else
            granFactor = data.data(:,140);
        end
    elseif polychoric == 1
        data = importdata('39term_polychoric_static.xlsx');
        BAI = data.data(:,21);
        BDI = data.data(:,22);
        granFactor = data.data(:,71);
    else
        data = importdata('39term_Pearson_static.xlsx');
        BAI = data.data(:,21);
        BDI = data.data(:,22);
        granFactor = data.data(:,70);
    end
end

%% create scatter plots of significant relationships
if dataSet == 18
    if windowed == 1
        if SDs == 1
            figure;
            scatter1 = scatter(granFactor,mPHQ8,[],rgb('MediumPurple'),'filled');
            xlabel('granularity factor (SD estimates)');
            ylim([0 24]);
            xlim([-1.5 2]);
            yticks([4 8 12 16 20 24]);
            set(gca,'fontsize',fontsize)
            ylabel('mean self-reported depression (PHQ-8)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_mPHQ8_scatter_plot_window_SDs'];
            saveas(scatter1,filename,'tiff');
            
            figure;
            scatter2 = scatter(granFactor,mGAD7,[],rgb('Tomato'),'filled');
            xlabel('granularity factor (SD estimates)');
            ylim([0 21]);
            xlim([-1.5 2]);
            set(gca,'fontsize',fontsize)
            yticks([3 6 9 12 15 18 21]);
            ylabel('mean self-reported anxiety (GAD-7)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_mGAD7_scatter_plot_window_SDs'];
            saveas(scatter2,filename,'tiff');
        else
            figure;
            scatter1 = scatter(granFactor,mPHQ8,[],rgb('MediumPurple'),'filled');
            xlabel('granularity factor #2 (mean estimates)');
            ylim([0 24]);
            set(gca,'fontsize',fontsize)
            yticks([4 8 12 16 20 24]);
            ylabel('mean self-reported depression (PHQ-8)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_mPHQ8_scatter_plot_window_Ms'];
            saveas(scatter1,filename,'tiff');
            
            figure;
            scatter2 = scatter(granFactor,mGAD7,[],rgb('Tomato'),'filled');
            xlabel('granularity factor #2 (mean estimates)');
            ylim([0 21]);
            set(gca,'fontsize',fontsize)
            yticks([3 6 9 12 15 18 21]);
            ylabel('mean self-reported anxiety (GAD-7)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_mGAD7_scatter_plot_window_Ms'];
            saveas(scatter2,filename,'tiff');
        end
    elseif polychoric == 1
        figure;
        scatter1 = scatter(granFactor,mPHQ8,[],rgb('MediumPurple'),'filled');
        xlabel('granularity factor #2 (overall estimates)');
        ylim([0 24]);
        set(gca,'fontsize',fontsize)
        yticks([4 8 12 16 20 24]);
        ylabel('mean self-reported depression (PHQ-8)');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_mPHQ8_scatter_plot_polychoric'];
        saveas(scatter1,filename,'tiff');

        figure;
        scatter2 = scatter(granFactor,mGAD7,[],rgb('Tomato'),'filled');
        xlabel('granularity factor #2 (overall estimates)');
        ylim([0 21]);
        set(gca,'fontsize',fontsize)
        yticks([3 6 9 12 15 18 21]);
        ylabel('mean self-reported anxiety (GAD-7)');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_mGAD7_scatter_plot_polychoric'];
        saveas(scatter2,filename,'tiff');
    else
        figure;
        scatter1 = scatter(granFactor,mPHQ8,[],rgb('MediumPurple'),'filled');
        xlabel('granularity factor #2 (overall estimates)');
        ylim([0 24]);
        set(gca,'fontsize',fontsize)
        yticks([4 8 12 16 20 24]);
        ylabel('mean self-reported depression (PHQ-8)');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_mPHQ8_scatter_plot'];
        saveas(scatter1,filename,'tiff');
    end
elseif dataSet == 39
    if windowed == 1
        if SDs == 1
            figure;
            scatter1 = scatter(granFactor,BAI,[],rgb('SeaGreen'),'filled');
            xlabel('granularity factor (SD estimates)');
            ylim([0 63]);
            set(gca,'fontsize',fontsize)
            yticks([9 18 27 36 45 54 63]);
            ylabel('self-reported anxiety (BAI)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_BAI_scatter_plot_window_SDs'];
            saveas(scatter1,filename,'tiff');
        else
            figure;
            scatter1 = scatter(granFactor,BAI,[],rgb('SeaGreen'),'filled');
            xlabel('granularity factor #2 (mean estimates)');
            ylim([0 63]);
            set(gca,'fontsize',fontsize)
            yticks([9 18 27 36 45 54 63]);
            ylabel('self-reported anxiety (BAI)');
            h1 = lsline;
            h1.Color = rgb('Gray');
            filename = ['dataSet' num2str(dataSet) '_BAI_scatter_plot_window_Ms'];
            saveas(scatter1,filename,'tiff');
        end
    elseif polychoric == 1
        figure;
        scatter1 = scatter(granFactor,BAI,[],rgb('SeaGreen'),'filled');
        xlabel('granularity factor #1 (overall estimates)');
        ylim([0 63]);
        set(gca,'fontsize',fontsize)
        yticks([9 18 27 36 45 54 63]);
        ylabel('self-reported anxiety (BAI)');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_BAI_scatter_plot_polychoric'];
        saveas(scatter1,filename,'tiff');
    else
        figure;
        scatter1 = scatter(granFactor,BAI,[],rgb('SeaGreen'),'filled');
        xlabel('granularity factor #1 (overall estimates)');
        ylim([0 63]);
        set(gca,'fontsize',fontsize)
        yticks([9 18 27 36 45 54 63]);
        ylabel('self-reported anxiety (BAI)');
        h1 = lsline;
        h1.Color = rgb('Gray');
        filename = ['dataSet' num2str(dataSet) '_BAI_scatter_plot'];
        saveas(scatter1,filename,'tiff');
    end
end

close all;
