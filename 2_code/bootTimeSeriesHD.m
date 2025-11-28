function dataBSHD = bootTimeSeriesHD(tp, exp, nsample)

% Bootstrap a time-series dataset. Closer timepoints have higher weights.

% tp            time points
% exp           expression data for all genes at each time point
% nsample       how many point to sample. Default: 100. 

% Jessline Haruman. Last update: 2025-08-06

%linspace on time
BS = linspace(0, tp(end), nsample).';

%set array for bootstraped gene expression
nGenes = width(exp);
BSExp = zeros(length(BS), nGenes);


%loop over genes
for g = 1:nGenes
  
    %loop over BS
    for i = 1:length(BS)
        
        %if time bootstraped is the match original data time, get expression value
        if ismember(BS(i), tp)
            inx = find(tp == BS(i), 1);  %find the index
            BSExp(i, g) = exp(inx, g);   
            
          
        %if not match, create a new value for gene expression
        else
       
            %calculate the distance between tp, closer point get 2.5 weight
            dist = 1./ (abs(BS(i)-tp).^2.5); 
    
            %normalization
            weight = normalize (dist, 'norm', 1); 
            bsd = randsample(exp(:,g), 1000, true, weight); 
            bsd = sum (bsd)/1000; 
    
            %assign new value
            BSExp (i, g) = bsd;
        end

     end
 
 
dataBSHD = [BS, BSExp];


end
