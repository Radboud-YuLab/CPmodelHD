function runRandG(numRuns, lambdaList, epsilonList, numGenes, log2_fc, gene_names, tp, fileName)
% Run a number of getCPsHD with different number of genes and genes. The result of
% this script will be used to determine the best lambda and epsilon.
% 
% input:
% numRums               Number of runs
% lambdaList            List consist of lambda values to test    
% epsilonList           List of small adjustment factor for non-negative linear regression
% numGenes              Number of genes to test
% log2_fc               Log 2 fold gene expression value
% gene_names            Gene names
% tp                    Experiment time point
% lambdaResultPath      Path to store all_result file in .txt or .csv
%
% Jessline H. Last update: 2025-11-30.

%%

% for saving result
fid = fopen(fileName, 'w');
fprintf(fid, 'numRuns\tgenes\tlambda\tepsilon\tnumCPs\tfinalCPs\tflag\trss\tprss\tr2\n');
fclose(fid);

% start loop
for r = 1:numRuns

    % random pick N genes
    randIndices = randperm(size(log2_fc, 1), numGenes);
    randGenesNames = gene_names(randIndices, :);
    
    genesStr = strjoin(randGenesNames, ';');

    
    % get corresponding data indices from the dataset
    exp = log2_fc(randIndices, :)';
    data = cat (2, tp', exp);
    
    % bootstrap on time
    dataBSHD = bootTimeSeriesHD(tp', exp, 100);
    
    % status 
    fprintf('Run %d / %d with genes: %s\n', r, numRuns, genesStr);

    
    for k = 1:length(lambdaList)
        lambda = lambdaList(k);
    
        % get tau for data bootstraping and rescaling
        % set interval for rescaling
        scale_intv = [0, max(tp)].* lambda; 
        [tauHD, ~] = getTauHD(dataBSHD, true, scale_intv);
        
        % set maxnCP
        % maxnCP = floor((length(tp)-2)/3)+2;
        maxnCP = 6; % user input, change if needed
    
        for e = 1:length(epsilonList)
            epsilon = epsilonList(e);
    
     
            % getCP
            output_cp_hd = getCPsHD(dataBSHD, maxnCP, tauHD, epsilon, data);
          
            % save output to csv 
            % extract row for saving the data to the csv
            
            numRows = size(output_cp_hd, 1);
        
            for i = 1:numRows
        
                numCPs = output_cp_hd{i,1};
                flag = output_cp_hd{i,3};
                rss = num2str(round(output_cp_hd{i,4},2));
                prss = num2str(round(output_cp_hd{i,5},2)); 
                r2 = num2str(round(output_cp_hd{i,6},2)); 
          
                CPCoords = output_cp_hd{i, 2};
        
                CPStr = strings(size(CPCoords, 1), 1);
                
                for j = 1:size(CPCoords,1)
                    CPStr(j) = strjoin(string(CPCoords(j,:)), ';');  % store each row as a cell entry
                end
                
                finalCPStr = strjoin(CPStr, ' | ');
            
                
                % append result to csv
                fid = fopen(fileName, 'a');
                fprintf(fid, '%d\t%s\t%.8f\t%.8f\t%d\t%s\t%d\t%s\t%s\t%s\n', ...
                    r, genesStr, lambda, epsilon, numCPs, finalCPStr, flag, rss, prss, r2);
                fclose(fid);
             end 
        end
    end
end
