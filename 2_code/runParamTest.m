% Main script to find best parameters values

% Jessline H. Last update: 2025-11-30
%%
numRuns = 20; %increasing the number of run might remove the inconsistency of optimal lambda !!! not always works ???

%linCutoff = 0.95;
lambdaList = [0.25, 0.16, 0.1, 0.5, 1, 2, 4 ,6];

%R2Cutoff = 0.6;
epsilon = [0.1, 0.2];  %adjustment factor, arbitary number between 0-1

% number random genes to select
numGenes = 2;
%% 
% find optimal parameters
lambdaResultFile = '4_processed_data\pom_main\resultHD\bestParam_1.csv';

% optional: graph generation
graphPath = '6_results\figuresCP_HD\miniTestLambda4';

runRandG(numRuns, lambdaList, epsilon, numGenes, pom_log2_fc, pom_gene_names, pom_days, ...
    lambdaResultFile);

[bestResult, all_results] = getBestParameters(lambdaResultFile);


%%
% optional: boxplot generation
% adjust the y axis if needed 

% input file
boxPlotPathInput    = '4_processed_data\pom_main\resultHD\bestParam_1.csv';
% output folder for image generated
boxPlotPathOutput   = '4_processed_data\pom_main\resultHD\boxplot\uhlitz_res14';

drawBoxPlot(boxPlotPathInput, boxPlotPathOutput);

