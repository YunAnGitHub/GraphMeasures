function [ local_efficiency ] = measure_local_efficiency( gPPI_weighted_matrix )
%MEASURE_LOCAL_EFFICIENCY Summary of this function goes here
%   Detailed explanation goes here


% 2018-06-28 Yun-An Huang
% the script is used to claculate the local efficiency by using the
% formula introduced by yu wang et al. 2017
% gPPI_weighted_matrix is a directed and weighted matrix.

% parameter

N_node = size(gPPI_weighted_matrix,1);
Max_w = max(gPPI_weighted_matrix(:));
Eff_array = zeros(1,N_node);

for itemp = 1:N_node
    
    numerator = 0;
    denominator = 0;
    
    % find neighbor of node i
    isnot_neighbor_idx = find( ~( gPPI_weighted_matrix(itemp,:) | gPPI_weighted_matrix(:,itemp).'));
    neighbor_idx = find( ( gPPI_weighted_matrix(itemp,:) | gPPI_weighted_matrix(:,itemp).'));
    
    for jtemp = neighbor_idx % for the neighbor
        
        for htemp = neighbor_idx % for the neighbor
            
            if jtemp ~= htemp                           
                
               gPPI_weighted_matrix_adapted = gPPI_weighted_matrix;
               gPPI_weighted_matrix_adapted(jtemp,htemp) = gPPI_weighted_matrix(itemp,jtemp)*gPPI_weighted_matrix(itemp,htemp)*gPPI_weighted_matrix(jtemp,htemp)/(Max_w^3);% replace weighted of j-> h
               gPPI_weighted_matrix_adapted(itemp,:) = 0; % remove node i, set node i as disconnected node
               gPPI_weighted_matrix_adapted(:,itemp) = 0; % remove node i, set node i as disconnected node
               gPPI_weighted_matrix_adapted(isnot_neighbor_idx,:) = 0; % remove the non-neightbor node, set node i as disconnected node
               gPPI_weighted_matrix_adapted(:,isnot_neighbor_idx) = 0; % remove the non-neightbor node, set node i as disconnected node
                                             
               [dist_array, track_cell] = dijkstra_weighted( 1./gPPI_weighted_matrix_adapted, jtemp );
                
               adapted_weighted = 1/dist_array(htemp); % the dist_array(htemp) is the shortest distance from jtemp to htemp in the adapted gPPI weighted matrix.
               
               numerator = numerator+ gPPI_weighted_matrix(itemp,jtemp)*gPPI_weighted_matrix(itemp,htemp)*adapted_weighted;
               denominator = denominator + gPPI_weighted_matrix(itemp,jtemp)*gPPI_weighted_matrix(itemp,htemp);
                
            end
            
            
        end
    end
    
    Eff_array(itemp) = numerator/denominator;
    
end

local_efficiency = sum(Eff_array)/N_node;


end

