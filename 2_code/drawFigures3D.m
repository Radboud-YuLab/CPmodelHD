% Generate a 3D plot for CP HD output. 

% This script includes a loop to automate the generation of 3D graph for 
% all CP output and it does not set a specifiv view for the graph. 
% Users can remove the loop to create a separate figures for each CP and 
% adjust the view point. 

% Jessline Haruman. Last update: 2025-11-29.

% set graphPath
graphPath = '6_results\figuresCP_HD\20251003';

% make directory
currentPath = pwd;
if ~isfolder(graphPath)
    mkdir(graphPath)
end

cd (graphPath)


%%
%draw & save graphs
for ngraphs = 1:height(output_cp_hd)
    P = graph_points(output_cp_hd{ngraphs, 1}, output_cp_hd{ngraphs, 2});
    gcf = figure('visible','off');
   

    % plot graph points
    plot3(P(:,1), P(:,2), P(:,3));
        
    hold on
     
    % plot data 
    scatter3(data(:,1), data(:,2), data(:,3), 10, "blue", "filled", "o");

    % plot bootstrapped data
    scatter3(dataBSHD(:,1), dataBSHD(:,2), dataBSHD(:,3), 10, "black", "o");

    % plot CPs
    P_CP = output_cp_hd{ngraphs, 2};
    scatter3(P_CP(:,1), P_CP(:,2), P_CP(:,3), 15, "red", "filled", "o");
 
    grid on
    grid minor
    
    ax = gca;

    maxOutputCPAxis = max(CPs);
    minOutputCPAxis = min(CPs);
    
    gtitle = ['dummydata2', '_', num2str(output_cp_hd{ngraphs, 1})];
        
    % legend
    title(gtitle, 'Interpreter', 'none')
    xlabel ('Time')
    ylabel ('Gene A')
    zlabel ('Gene B')
    hold off
    
    % save figure
    saveas(gcf,gtitle,'png')
    
end