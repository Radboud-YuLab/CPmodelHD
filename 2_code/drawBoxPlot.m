function [] = drawBoxPlot(boxPlotPathInput, boxPlotPathOutput)
% Generate boxplot to visualize the median results for X times algorithm runs, showing the distribution of 
% each goodness of fit parameters median
%
% input:
% boxPlotPathInput          Path to input file in .txt or .csv containing the median values for all the 
%                           goodnes of fit parameters along with their corresponding lambda values
%
% output:
% boxPlotPathOutput          Path to store generated boxplot
%
% Jessline H. Last update: 2025-11-30

% load file
data = readtable(boxPlotPathInput, 'FileType','text','Delimiter','\t');

% check output folder
% make directory if folder not exist

if ~isfolder(boxPlotPathOutput)
    mkdir(boxPlotPathOutput)
end

% get number of CPs
cpList = unique(data.numCPs);

% metrics to plot
metrics = {'rss', 'prss', 'r2'};

%%
for i = 1:length(cpList)
    
    cp = cpList(i);
    
    %filter data based on number of CP
    dataSub = data(data.numCPs == cp, :);

    %extract column
    lambda = dataSub.lambda;
    [uniqueLambdas, ~, lambdaGroup] = unique(lambda);

    %convert to categorial
    lambdaCat = categorical(lambda, uniqueLambdas);
    
    for m = 1:length(metrics)
        currMetrics = metrics{m};
        metrics_data = dataSub.(currMetrics);
        
        %plot
        figure;
        hold on;
        
        boxchart(lambdaCat, metrics_data, 'BoxFaceColor', [0 0.4470 0.7410]);
        
        %overlay datapoints
        for j = 1:length(uniqueLambdas)
        
        x = j + 0.15 * (rand(sum(lambdaGroup == j), 1) - 0.5); % jitter
        y = metrics_data(lambdaGroup == j);
        scatter(x, y, 20, 'MarkerEdgeColor', 'r');

        end

        % plot legend
        xlabel('Lambda');
        ylabel(currMetrics);
        title(sprintf('%s distribution (nCPs = %d)', currMetrics, cp));
        

        % set y-limits (customize per metrics if needed)
        
        if ismember({'r2'}, currMetrics)
            ylim([0 1]);

        else
            ylim([0 200]);
        end
                
        grid on;
        hold off;
        
        %save figures
        currTime = datetime('now', 'Format', 'yyyyMMdd_HHmmSS');
        fileName = fullfile(boxPlotPathOutput, sprintf('%s_%s_boxplot_CP%d.png', currTime, currMetrics, cp) );

        saveas(gcf, fileName);
        close(gcf);

    end
end
end
