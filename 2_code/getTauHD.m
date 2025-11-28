function [tauHD, scale_data] = getTauHD(data, rescale_data, lambda)

% Calculate tau in high dimension for CP model fitting. 

% data              An n x m array. First row is the time of experiments, 
%                   other rows are the genes. Every column is gene
%                   expression data corresponding to the timepoints.
% rescale_data      Boolean. To rescale data or not.
% lambda            Scaling factor

% Jessline Haruman. Last update: 2025-11-28


%rescaling 

scale_data = data;

if rescale_data 
    % rescale gene based on time value
    tp_min = min(scale_data(:, 1));
    tp_max = max(scale_data(:, 1));

    for col_gene = 2:size(scale_data, 2) % time column skipped
        scale_intv = [0, (tp_max - tp_min)] .* lambda;
        
        scale_data(:, col_gene) = rescale(scale_data(:, col_gene), scale_intv(1), scale_intv(2));
    end
else
    disp("data is not rescaled")

end
   


% calculate distance

distHD = sqrt(sum(diff(scale_data(:, 1:end)).^2, 2));
 
% normalize distance

norm_distHD = cumsum(distHD)/sum(distHD);

% append tau

tauHD = cat (1, 1e-5, norm_distHD(1:length(norm_distHD)-1,:), (1 - 1e-5));






