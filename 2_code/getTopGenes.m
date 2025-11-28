function topGenes = getTopGenes(output_cp_hd, geneNames, numGenes)

% Get genes with elevated expression values in each CP time point exclude
% for the first and the last data points

% output_cp_hd          An array of n x 6 array. Second column is important
%                       for this function (CP coordinate)
% geneNames             An array of 1 x m gene names
% numGenes              Number of top genes. Default: 10

% Jessline Haruman. Last update: 2025-08-06

if nargin < 3 
    numGenes = 10;
end

% transpose the gene name into a column vector if they are in a row format
if size(geneNames,2) > 1 
    geneNames = geneNames.'; 
end

% output structure
topGenes = {}; 

for n = 1:size(output_cp_hd,1)
    numCP = output_cp_hd{n,1};

    if numCP == 2
        continue;
    end

    CPCoord = output_cp_hd{n,2};
    
    % error handling: if no CP value is empty
    if isempty(CPCoord)
            continue;
    end

    CPtime = CPCoord(:,1);
    CPexp = CPCoord(:,2:end);
    
    
    nCPtime = size(CPtime,1); 

    for i = 2:nCPtime-1 % exclude first and last datapoint
        values = CPexp(i,:);
        [~, idx] = sort(values, 'descend');
        topGenesNames = geneNames(idx(1:numGenes));
        geneList = strjoin(topGenesNames, ',');

        topGenes(end+1,:) = {numCP, CPtime(i), geneList};

    end

end



