function output_cp_hd = getCPsHD(dataBSHD, maxnCP, tauHD, epsilon, data)

% Fitting CP models to the given data.

% dataBSHD          Bootstraped data.
% maxnCP            Maximum number of CP.
% tauHD             An array m x 1 dimension. 
% epsilon           Small arbitary number to prevent overlapping between CP
% data              An n x m array. First row is the time of experiments, 
%                   other rows are the genes. Every column is gene.
%                   expression data corresponding to the timepoints.

% output_cp         An array of class cell. 
%                   First column is the number of CPs.
%                   Second columns is array of CPs coordinate (x,y). 
%                   Each row will contain number of CPs and their coordinate accordingly.

% Jessline Haruman. Last update: 2025-07-23

% set up output structures
CPs = [];
output_cp_hd (1,1) = num2cell (2);
output_cp_hd (1,2) = {CPs};

% loop over maxnCP
for m = 1:(maxnCP-1)

    % reset cp for each itteration
    CPs = [];


    % for maxnCP=2, use first and last data points
    if m == 1
        CPs (1,:) = dataBSHD(1,:);
        CPs (2,:) = dataBSHD(end,:);
   
    else
        
        % calculate bernstein basis for nCPs
        B_all = bernstein_basis (m, tauHD); % bernstein_basis(n, tau) // n = degree of bernstein basis (number of control points - 1)
        
        % remove first and last basis in B
        G = B_all;
        G(:,1) = []; %first basis 
        G(:,end) = []; %last basis 
        
        
        % find CP for time with no negative value

        % exclude the last datapoint from fitting 
        data_t = dataBSHD(:, 1) - (B_all(:,end).*dataBSHD(end,1));    
        
        % calculate basis for non-negative linear regression
        [w, n] = size(G);           %w = n rows, c = n columns
        B_t = zeros(w,n);

        for i = 1:n 
            B_t(:, i) = sum(G(:, i:n), 2);
        end 
        
   
        % non-negative linear regression
        % set lower bound, enforce D > 0 
        lb = epsilon * ones(n, 1) * max(data(:,1)); %tp 

        D = lsqlin(B_t, data_t, [], [], [], [], lb, []);  % Solves min ||Cx - d||²

        % get time CP value from D
        CP_time = zeros(n, 1);       
        CP_time (1) = D(1);

        for i = 2:n
            CP_time(i) = CP_time(i-1) + D(i);
        end
       
       
        % multilinear regression 

        % find CP for genes
        for i = 2:size(dataBSHD,2) 
            data_n = dataBSHD(:, i) - ((B_all(:,1)).*dataBSHD(1,i) + B_all(:,end).*dataBSHD(end,i));
           
            [bn, ~, ~] = regress (data_n, G);
            CPs = cat(2, CPs, bn);
        end
        
        % add CP time 
        CPs = cat(2, CP_time, CPs);
        
        % add CP0 and CPn from the first and last values of data         
        CPs = cat (1, dataBSHD(1,:), CPs, dataBSHD(end,:)) ;
    
    end 
    
    % append CPs to output
    
output_cp_hd(m,1) = num2cell (m+1); %number of CP
    output_cp_hd(m,2) = {CPs};
end

%%

% analyzing the CP output
for j = 1:height(output_cp_hd)
    nCP = cell2mat(output_cp_hd(j,1));
    b = output_cp_hd{j,2}(:,1); %time column
    
    %flag 
    %0 nothing wrong with it 
    %1 not monotonously increasing
    %2 time is negative
    %3 both problem

    % both problem
    if ~all(diff(b)>0) && any(b < 0)
        flag = 3;

    % if CPs on x (time) are not monotonously increasing: 
    elseif ~all(diff(b)>0)
        flag = 1;
    
    % if time is negative
    elseif any(b < 0) 
        flag = 2;
    
    % no flag
    else 
        flag = 0;

    end
    
    output_cp_hd(j,3) = {flag};
    
    
    % convert Bernstein basis control points to power basis coefficients.
    n = changeBasis(b);
    
    % an empty array for predicted y-values
    % -2 to remove first and last data points, will be added back later
    
    nGene = width(data)-1;
    nDataPoints = height(data);

    y_hat = zeros((nDataPoints-2),(nGene));

    for k = 2:(nDataPoints-1)
        % calculate r (tau) from x-value (time)
        % first and last data points are not counted here (tau = 0 & tau = 1) 
        poly = n;
        poly(end) = poly(end) - data(k,1);
        r = roots(poly);
        r = real(r(imag(r) == 0));
        r = r(round(r,4)>=0 & round(r,4) <=1);
   
        for d = 1:(nGene) % exclude time column
            % use r to calculate fitted y-value (y_hat)
            P = calculate_y(nCP, output_cp_hd{j,2}(:,d+1), r);
            y_hat(k, d) = P;
        end
    end
    
    % add back first and last data points 
    y_hat = cat (1, y_hat, data(end,2:end));

    % calculate goodness of fit parameter
    %1. calculate RSS
    rss = sum(sum((data(:, 2:end) - y_hat).^2));
    
    %2. calculate nCP-penalized RSS
    sp_rss = rss * (nCP-1)^2;
    
    %3. calculate R squared
    sst = sum(sum((data(:,2:end) - mean(data(:,2:end), 1)).^2));
    rsq = 1-(rss/sst);
    if rsq < 0
        rsq = 0;
    end
    
    output_cp_hd(j,4) = num2cell(rss);
    output_cp_hd(j,5) = num2cell(sp_rss);
    output_cp_hd(j,6) = num2cell(rsq);
    
end 

end


% local
function B = bernstein_basis(n, tau)
    % Compute Bernstein basis polynomials up to degree n at point tau
    % n = number of control points - 1
    
    % Initialize matrix to store the basis polynomials
    B = zeros(length (tau), n+1);
    
    % Calculate the basis polynomials
    for i = 0:n
        B(:, i+1) = nchoosek(n, i) * tau.^i .* (1 - tau).^(n - i);
    end
end
