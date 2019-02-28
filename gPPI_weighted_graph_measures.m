function [ global_measures ] = gPPI_weighted_graph_measures( gPPI_weighted_matrix_ori )
%GPPI_GRAPH_MEASURES Summary of this function goes here
%   Detailed explanation goes here


% 2019-Feb-28 Yun-An Huang


% 2018-Jun-28 Yun-An Huang
% using local efficiency from yu wang et al., 2017.
% using cluster coefficient from Miyajima and Sakuragawa et el., 2014. 

% 2018-Jun-15 Yun-An Huang
% using another definition of local efficiency. (the previous version use
% the formula from Rubinov et al., 2010)

% 2018-Apr-17 Yun-An Huang
% set the diagonal as zero.

% Yun-An Huang 2018-March-21
% this function was adopted to use dijkstra_weighted_multi which report the
% multi tracks. This is important to calculate betweenness centrality.

% Yun-An Huang 2018-Feb-06
% this function is calculating the global index of a weighted gPPI matrix.
% noteworthy, gPPI_matrix here is weighted and directed.
% the global measures refer to M. Rubinov, O. Sporns "Complex network
% measures of brain connectivity: Uses and interpretations" (2010)

% input: 
% gPPI_weighted_matrix: the n by n directed weight matrix.

% output several graphical measures:
%
% 1. sum_of_weight: the sum of all weight of gPPI matrix.
% 2. out_degree: the out degree of nodes.
% 3. in_degree: the in degree of nodes.
% 4. ShortPathLength_matrix: the matrix with shortest path length between each two nodes.
% 5. track_shortpath: the cell array store the tracks of shortest path
%                     length between each two nodes.
% 6. num_triangles : the arrary with number of triangles of nodes
% 7. char_path_length: the characteristic path length of the graph
% 8. global_eff: the global efficiency of the graph
% 9. cluster_coeff: the cluster coefficiency of the graph
% 10. transitivity: the transitivity of the graph
% 11. local_efficiency: the local efficiency
% 12. modularity: the modularity
% 13. closeness_centrality: the closeness centrality
% 14. betweenness_centrality: the betweenness centrality
% 15. within_module_out_degree: the between module degree.
% 16. within_module_in_degree: the between module degree.
% 17. mean_weight_neighbor_degree: the average neighbor degree.
% 18. assortativity: the assortativity
% 19. small_worldness: the small wordness.



%% initial parameter


node_num = size(gPPI_weighted_matrix_ori,1);


% roudn the matrix. to the 5 digit decimal point
gPPI_weighted_matrix_round = round(gPPI_weighted_matrix_ori * 10000)/10000;
gPPI_weighted_matrix = zeros(node_num,node_num);

idx = eye(node_num); % the diag index;
gPPI_weighted_matrix(~idx) = gPPI_weighted_matrix_round(~idx); % remove diagonal data

global_measures = [];
global_measures.sum_of_weight = [];
global_measures.out_degree = [];
global_measures.in_degree = [];
global_measures.ShortPathLength_matrix = [];
global_measures.track_shortpath = [];
global_measures.num_triangles = [];
global_measures.char_path_length = [];
global_measures.global_eff = [];
global_measures.cluster_coeff = [];
global_measures.transitivity = [];
global_measures.local_efficiency = [];
global_measures.modularity = [];
global_measures.closeness_centrality = [];
global_measures.betweenness_centrality = [];
global_measures.within_module_out_degree = [];
global_measures.within_module_in_degree = [];
global_measures.mean_weight_neighbor_degree = [];
global_measures.assortativity = [];
global_measures.small_worldness = [];

%% Sum of weights

global_measures.sum_of_weight = sum(gPPI_weighted_matrix(:));

%% out_degree

global_measures.out_degree = sum(gPPI_weighted_matrix,2)';

%% in_degree

global_measures.in_degree = sum(gPPI_weighted_matrix,1);

%% Shortest weighted path length, SPL 


distance_matrix = 1./gPPI_weighted_matrix;
global_measures.ShortPathLength_matrix = zeros(node_num,node_num); % the matrix store the path length of each two nodes.
global_measures.track_shortpath = [];

for itemp = 1:node_num
    
    [dist_array, track_cell] = dijkstra_weighted_multi(distance_matrix,itemp); % for each node becaome a source
    global_measures.ShortPathLength_matrix(itemp,:) = dist_array; 
    
    for jtemp = 1:node_num
        
        global_measures.track_shortpath{itemp,jtemp} = track_cell(jtemp,:); % the minimum path track from node i to node j
    end
    
    
end



%% number of triangles,  the array store the number of triangles around node i.

global_measures.num_triangles = measure_num_triangles(gPPI_weighted_matrix);


%% characteristic path length.


% the diagonal of ShortPathLength_matrix should be zero
global_measures.char_path_length = sum(global_measures.ShortPathLength_matrix(:))/node_num/(node_num-1); 


%% global efficiency
efficiency_matrix = 1./global_measures.ShortPathLength_matrix ;

for itemp = 1:node_num
    efficiency_matrix(itemp,itemp)=0; % do not consider the diagonal 
end

global_measures.global_eff = sum(efficiency_matrix(:))/node_num/(node_num-1);


%% clustering coefficient.

% using the hramonic mean introduced by miyajima and sakuragawa.

global_measures.cluster_coeff = measure_clustering_coefficient(gPPI_weighted_matrix);

%% Transitivity

trans_arr_temp = zeros(1,node_num);
out_degree = global_measures.out_degree;
in_degree = global_measures.in_degree;

for itemp = 1:node_num

    trans_arr_temp(itemp) = ( (out_degree(itemp)+in_degree(itemp))*(out_degree(itemp)+in_degree(itemp)-1)-2*(gPPI_weighted_matrix(itemp,:)*gPPI_weighted_matrix(:,itemp)));

end

global_measures.transitivity = sum( global_measures.num_triangles) / sum(trans_arr_temp);


%% local efficiency
% since this is a weighted matrix withiout setting a threshold, every node is the neighbor of others.

% using the local efficiency introduced by yu wang et al., 2017


global_measures.local_efficiency = measure_local_efficiency(gPPI_weighted_matrix);


%% modularity
% since this is a weighted matrix without setting a threshold, every nodes
% is in the same module.

modularity_matrix_temp = zeros(node_num,node_num);

for itemp = 1:node_num
    for jtemp = 1:node_num
    
        modularity_matrix_temp(itemp,jtemp) = gPPI_weighted_matrix(itemp,jtemp)-(global_measures.out_degree(itemp)*global_measures.in_degree(jtemp)/global_measures.sum_of_weight);
        
    end
end

global_measures.modularity = sum(modularity_matrix_temp(:))/global_measures.sum_of_weight;


%% closeness centrality

closeness_centrality = zeros(1,node_num);



for itemp = 1:node_num
   
    closeness_centrality(itemp) = (node_num-1)/sum(global_measures.ShortPathLength_matrix(itemp,:));    % the self connection of short path length is zero.
    
end

global_measures.closeness_centrality = closeness_centrality;

%% Betweeness centrality
% the betweeness centrality is the portion of shortest paths pass through
% node i in the all of the shortest paths.


betweenness_centrality = zeros(1,node_num);

for itemp = 1:node_num
    
    betweenness_centrality(itemp) = measure_betweeness_centrality(global_measures.track_shortpath,itemp); % calculate the betweenness centrality of node i.
    
end

global_measures.betweenness_centrality = betweenness_centrality;
%% Within-module degree z-score.
% since we have directed weighted graph, such that all nodes are in the same module. 
% the degree separate into out-degree and in-degree. 
% 

global_measures.within_module_out_degree = (global_measures.out_degree-mean(global_measures.out_degree(:)))/std(global_measures.out_degree(:));

global_measures.within_module_in_degree = (global_measures.in_degree-mean(global_measures.in_degree(:)))/std(global_measures.in_degree(:));

%% Particiation coefficient
% since there is only one module, the participation coefficient is 0.

%% network motifs
% since we used the weighted gPPI matrix without threshold, every links are
% existed, it is meaningless to calculate motif since the network motif would be the same as
% in a random graph.

%% degree distribution.
% i am not sure what is the definition of probability of degree k'
% is it coming from the random network?

%% Average neighbor degree
%

mean_weight_neighbor_degree = zeros(1,node_num);

for itemp = 1: node_num

    sum_neighbor = 0;
    for jtemp = 1:node_num
        
        sum_neighbor = sum_neighbor+(gPPI_weighted_matrix(itemp,jtemp)+gPPI_weighted_matrix(jtemp,itemp))*(global_measures.out_degree(jtemp)+global_measures.in_degree(jtemp));

    end
    
    mean_weight_neighbor_degree(itemp) = sum_neighbor/(2*global_measures.out_degree(itemp)+global_measures.in_degree(itemp));
end

global_measures.mean_weight_neighbor_degree = mean_weight_neighbor_degree;

%% assortativity coefficient
%

l_inv =  1/global_measures.sum_of_weight;

wkij_multiply = 0;

wkij_sum = 0;

wkij_sum_square = 0;

for itemp = 1: node_num
    for jtemp = 1: node_num
        
%         wkij_multiply = wkij_multiply+  global_measures.out_degree(itemp) * global_measures.in_degree(jtemp);
%         
%         wkij_sum = wkij_sum+ 0.5 *(global_measures.out_degree(itemp) + global_measures.in_degree(jtemp));
%     
%         wkij_sum_square = wkij_sum_square + 0.5 * (global_measures.out_degree(itemp)^2 + global_measures.in_degree(jtemp)^2);

   
        wkij_multiply = wkij_multiply+ gPPI_weighted_matrix(itemp,jtemp)* global_measures.out_degree(itemp) * global_measures.in_degree(jtemp);
        
        wkij_sum = wkij_sum+ 0.5 * gPPI_weighted_matrix(itemp,jtemp)*(global_measures.out_degree(itemp) + global_measures.in_degree(jtemp));
    
        wkij_sum_square = wkij_sum_square + 0.5 * gPPI_weighted_matrix(itemp,jtemp)*(global_measures.out_degree(itemp)^2 + global_measures.in_degree(jtemp)^2);
    end
end

k_term_1 = l_inv * wkij_multiply ;

k_term_2 = ( l_inv * wkij_sum )^2;

k_term_3 = l_inv * wkij_sum_square;

global_measures.assortativity = (k_term_1 - k_term_2)/( k_term_3 - k_term_2);

%% Small worldness

global_measures.small_worldness = global_measures.cluster_coeff  / global_measures.char_path_length ; 





