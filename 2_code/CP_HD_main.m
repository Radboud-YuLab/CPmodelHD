% Main script for data analysis CP HD

% Jessline Haruman. Last update: 2025-11-29

%%
close all
clear

%% Pommerenke 2012 dataset, viral infection-induced immune response
% read in the Pommerenke 2012 dataset
pom_days = [0, 1, 2, 3, 5, 8, 10, 14, 18, 22, 26, 30, 40, 60];
[pom_gene_names, pom_log2_fc, pom_clusters] = readPommerenke;

%%
% tune/set data-specific parameters
lambda = 0.005;
epsilon = 0.001;

tp = pom_days.';
exp = pom_log2_fc.';
data = cat (2, tp, exp);

% fit CP model for all genes

% bootstrap on time
dataBSHD = bootTimeSeriesHD(tp, exp, 100);

% get tau for data bootstraping and rescaling
% set interval for rescaling
[tauHD, ~] = getTauHD(dataBSHD, true, lambda);
   
% set maxnCP 
maxnCP = floor((length(tp)-2)/3)+2;

% get CP HD
output_cp_hd = getCPsHD(dataBSHD, maxnCP, tauHD, epsilon, data);


%% Uhlitz 2017 dataset, ERK signaling
% read in Uhlitz 2017 dataset
uhl_hrs = [0,0.5,1,2,3,4,6,8,10];
[uhl_gene_names, ~, uhl_log2_fc] = readUhlitzTimeSeries;

%%
% tune/set data-specific parameters
lambda = 0.2;
epsilon = 0.0001;

tp = uhl_hrs.';
exp = uhl_log2_fc.';
data = cat (2, tp, exp);

% fit CP model for all genes

% bootstrap on time
dataBSHD = bootTimeSeriesHD(tp, exp, 100);

% get tau for data bootstraping and rescaling
% set interval for rescaling
[tauHD, ~] = getTauHD(dataBSHD, true, lambda);
   
% set maxnCP 
%maxnCP = floor((length(tp)-2)/3)+2;
maxnCP = 6;

% get CP HD
output_cp_hd = getCPsHD(dataBSHD, maxnCP, tauHD, epsilon, data);

%%
% validation
topGenes = getTopGenes(output_cp_hd, uhl_gene_names, 10);

 

%% Qu 2018 dataset, keratinocyte differentiation
% read in the Qu 2018 keratinocyte differentiation dataset
 
qu_days = [0,1,2,3,4,5,6,7];
[qu_genes, expr_wt, ~] = readQu;

% tune/set data-specific parameters
lambda = 0.25;
epsilon = 0.01;

tp = qu_days.';
exp = expr_wt.';
data = cat (2, tp, exp);

% fit CP model for all genes
tic
% bootstrap on time
dataBSHD = bootTimeSeriesHD(tp, exp, 100);
toc
% get tau for data bootstraping and rescaling
% set interval for rescaling
scale_intv = [0, max(tp)].* lambda; 
[tauHD, scale_data] = getTauHD(dataBSHD, true, scale_intv);
   
% set maxnCP 
maxnCP = floor((length(tp)-2)/3)+2;

% get CP HD
output_cp_hd = getCPsHD(dataBSHD, maxnCP, tauHD, epsilon, data);
