function [bestResult, all_results] = getBestParameters(fileName,graphPath)

% Determine the optimal lambda value by running the algorithm X number of times. For each run, it calculates the goodness
% of fit parameters to get best lambda and epsilon

% fileName          File containing the CP coordinates and the goodness of
%                   fit parameters for several number of runs 

% Jessline H. Last update: 2025-11-30.

% load files
data = readtable(fileName, 'FileType','text','Delimiter','\t'); 

% get the top 2 highest number of CP 
uniqueNumCPs = sort(unique(data.numCPs), 'descend');
nTop = min(3, numel(uniqueNumCPs));
topNumCPs = uniqueNumCPs(1:nTop);


dataSub = data(ismember(data.numCPs, topNumCPs), :);

% get unique pairs
pairs = unique(round([dataSub.lambda, dataSub.epsilon],12), 'rows', 'stable');
nPairs = size(pairs,1);

% result table ouput
all_results = table( ...
    pairs(:,1), pairs(:,2), ...
    nan(nPairs,1), nan(nPairs,1), nan(nPairs,1), nan(nPairs,1), ...
    'VariableNames', {'lambda','epsilon','median_R2','median_pRSS','median_RSS','score'});

% get median
for i = 1:nPairs
    lam = pairs(i,1);
    eps = pairs(i,2);

    rows = dataSub.lambda == lam & dataSub.epsilon == eps;

    all_results.median_R2(i)   = median(dataSub.r2(rows),   'omitnan');
    all_results.median_pRSS(i) = median(dataSub.prss(rows), 'omitnan');
    all_results.median_RSS(i)  = median(dataSub.rss(rows),  'omitnan');
end

% rank the median by sorting them
[~, rank_r2] = sort(all_results.median_R2, 'descend');
[~, rank_prss] = sort(all_results.median_pRSS, 'ascend');
[~, rank_rss] = sort(all_results.median_RSS, 'ascend');


% score values are based on the addition of index(rank) for each lambda for each median
% e.g for lambda 0.25, 
% r2 rank 1 is on the idx 3, 
% prss rank 1 is on the idx 1,
% and rss rank 1 is on the idx 2 
% then the total score would be 6. 
% after getting all the score for each lambda, then top 2 score to be selected
% the lower the score the better


for i = 1:nPairs
    all_results.score(i) = ...
        find(rank_r2 == i) + ...
        find(rank_prss == i) + ...
        find(rank_rss == i);
end

% get score
all_results.score = zeros(height(all_results),1);

% compute score: sum of ranks
for i = 1:height(all_results)
    idx_r2   = find(rank_r2 == i);
    idx_prss = find(rank_prss == i);
    idx_rss = find(rank_rss == i);

    all_results.score(i) = idx_r2 + idx_prss + idx_rss;
end

% sort by score
all_results = sortrows(all_results, 'score');

% best pair
bestResult = all_results{1, {'lambda','epsilon'}};




%%
% optional: generate 3D CP graph
% if not graph path provided then the algorithm will not generate graph

% graph can only be generated if the number of genes is 2 (3D graph)


if nargin >= 3 && ~isempty(graphPath) 
    % make directory if folder not exist
    
    if ~isfolder(graphPath)
        mkdir(graphPath)
    end

    cd (graphPath)


    for i = 1:height(dataSub)
        cp = cps{i,:};
        cpCoord = split(cp, "|");
    
        cp_matrix = zeros(length(cpCoord), 3);
        
        for k = 1:length(cpCoord)
            cp_matrix(k, :) = str2double(strsplit(cpCoord{k}, ';'));
        end
    
        P = graph_points(numcps(i), cp_matrix);
    
        gcf = figure('visible','off');
        
    
        % plot graph points
        plot3(P(:,1), P(:,2), P(:,3));
            
        hold on
    
        gene1 = dataSub.Gene1{i};
        gene2 = dataSub.Gene2{i};
    
        gene1_idx = find(ismember(gene_names, gene1));
        gene2_idx = find(ismember(gene_names, gene2));
        
        geneidx = [gene1_idx, gene2_idx];
        
        exp = log2_fc(geneidx, :)';
    
        dataPlot = cat (2, tp', exp);
    
        % plot data 
        scatter3(dataPlot(:,1), dataPlot(:,2), dataPlot(:,3), 10, "blue", "filled", "o");
        
    
        % plot CPs
        P_CP = cp_matrix;
        scatter3(P_CP(:,1), P_CP(:,2), P_CP(:,3), 15, "red", "filled", "o");
     
        grid on
        grid minor
           
        lambda = strrep(num2str(lambdas(i)), '.', '_');
    
        gtitle = [gene1, '_', gene2, '_n_', num2str(numcps(i)), '_lambda_', lambda];
        
        title(gtitle, 'Interpreter', 'none')
        xlabel ('Time')
        ylabel (gene1)
        zlabel (gene2)
        
        annotation('textbox',[.9 .5 .1 .2], ...
        'String', sprintf('r2 = %.2f', r2s(i)), 'EdgeColor','none')
    
        saveas(gcf,gtitle,'png')
    
    end

end 

end
%----------------------------------------------------------------------------------------------------------------------

% local
function P = graph_points(n, cp)
    % generate points from CPs for graphing
    % n = number of CPs
    
    n1 = n-1;

    % calculate the binomial coefficient (x!/(y!(x-y)!)) 
    for i = 0:1:n1
        binomial_coefficient(i+1) = factorial(n1)/(factorial(i)*factorial(n1-i));   
    end

    basis = [];
    UB = [];
    for tau = 0:0.01:1
        for k = 1:n
            UB(k) = binomial_coefficient(k)*((1-tau)^(n-k))*(tau^(k-1));
        end
        basis = cat(1,basis,UB);    
    end

    P = basis*cp;

end