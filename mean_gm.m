function [ global_measures ] = mean_gm( gm_cell )
%GPPI_GRAPH_MEASURES Summary of this function goes here
%   Detailed explanation goes here

% Yun-An Huang 2018-March-27
% update some of the global measures. 
% the measures of node scale are not included.

% Yun-An Huang 2018-Mar-02
% this function is calculating the mean global measures of gm.
% gm_cell is a cell;
%


total_num = length(gm_cell);

global_measures = [];

%% sum_of_weight

arr_sum_of_weight=zeros(1,total_num);

for p_temp = 1:total_num

    arr_sum_of_weight(p_temp) = gm_cell{p_temp}.sum_of_weight;

end

global_measures.sum_of_weight.mean = mean(arr_sum_of_weight(:));
global_measures.sum_of_weight.std = std(arr_sum_of_weight(:));
global_measures.sum_of_weight.all = arr_sum_of_weight(:);
%% char_path_length

arr_char_path_length=zeros(1,total_num);

for p_temp = 1:total_num

    arr_char_path_length(p_temp) = gm_cell{p_temp}.char_path_length;

end

global_measures.char_path_length.mean = mean(arr_char_path_length(:));
global_measures.char_path_length.std = std(arr_char_path_length(:));
global_measures.char_path_length.all = arr_char_path_length(:);

%% global_eff

arr_global_eff=zeros(1,total_num);

for p_temp = 1:total_num

    arr_global_eff(p_temp) = gm_cell{p_temp}.global_eff;

end

global_measures.global_eff.mean = mean(arr_global_eff(:));
global_measures.global_eff.std = std(arr_global_eff(:));
global_measures.global_eff.all = arr_global_eff(:);

%% cluster_coeff

arr_cluster_coeff=zeros(1,total_num);

for p_temp = 1:total_num

    arr_cluster_coeff(p_temp) = gm_cell{p_temp}.cluster_coeff;

end

global_measures.cluster_coeff.mean = mean(arr_cluster_coeff(:));
global_measures.cluster_coeff.std = std(arr_cluster_coeff(:));
global_measures.cluster_coeff.all = arr_cluster_coeff(:);

%% transitivity

arr_transitivity=zeros(1,total_num);

for p_temp = 1:total_num

    arr_transitivity(p_temp) = gm_cell{p_temp}.transitivity;

end

global_measures.transitivity.mean = mean(arr_transitivity(:));
global_measures.transitivity.std = std(arr_transitivity(:));
global_measures.transitivity.all = arr_transitivity(:);
%% local_efficiency

arr_local_efficiency=zeros(1,total_num);

for p_temp = 1:total_num

    arr_local_efficiency(p_temp) = gm_cell{p_temp}.local_efficiency;

end

global_measures.local_efficiency.mean = mean(arr_local_efficiency(:));
global_measures.local_efficiency.std = std(arr_local_efficiency(:));
global_measures.local_efficiency.all = arr_local_efficiency(:);


%% modularity

arr_modularity = zeros(1,total_num);

for p_temp = 1:total_num

    arr_modularity(p_temp) = gm_cell{p_temp}.modularity;

end

global_measures.modularity.mean =  mean( arr_modularity(:));
global_measures.modularity.std =  std( arr_modularity(:));
global_measures.modularity.all =   arr_modularity(:);

% %%  closeness centrality
% 
% arr_closeness_centrality = zeros(1,total_num);
% 
% for p_temp = 1:total_num
% 
%     arr_closeness_centrality(p_temp) = gm_cell{p_temp}.closeness_centrality;
% 
% end
% 
% global_measures.closeness_centrality =  mean( arr_closeness_centrality(:));
% 
% %%  Betweeness centrality
% 
% arr_betweenness_centrality = zeros(1,total_num);
% 
% for p_temp = 1:total_num
% 
%     arr_betweenness_centrality(p_temp) = gm_cell{p_temp}.betweenness_centrality;
% 
% end
% 
% global_measures.betweenness_centrality =  mean( arr_betweenness_centrality(:));
% 
% %%  Within-module outdegree
% 
% arr_within_module_out_degree = zeros(1,total_num);
% 
% for p_temp = 1:total_num
% 
%     arr_within_module_out_degree(p_temp) = gm_cell{p_temp}.within_module_out_degree;
% 
% end
% 
% global_measures.within_module_out_degree =  mean( arr_within_module_out_degree(:));
% 
% %%  Within-module indegree
% 
% arr_within_module_in_degree = zeros(1,total_num);
% 
% for p_temp = 1:total_num
% 
%     arr_within_module_in_degree(p_temp) = gm_cell{p_temp}.within_module_in_degree;
% 
% end
% 
% global_measures.within_module_in_degree =  mean( arr_within_module_in_degree(:));
% 
% %%  Average neighbor degree
% 
% arr_mean_weight_neighbor_degree = zeros(1,total_num);
% 
% for p_temp = 1:total_num
% 
%     arr_mean_weight_neighbor_degree(p_temp) = gm_cell{p_temp}.mean_weight_neighbor_degree;
% 
% end
% 
% global_measures.mean_weight_neighbor_degree =  mean( arr_mean_weight_neighbor_degree(:));

%% assortativity coefficient

arr_assortativity = zeros(1,total_num);

for p_temp = 1:total_num

    arr_assortativity(p_temp) = gm_cell{p_temp}.assortativity;

end

global_measures.assortativity.mean =  mean( arr_assortativity(:));
global_measures.assortativity.std =  std( arr_assortativity(:));
global_measures.assortativity.all =   arr_assortativity(:);

%% Small worldness

arr_small_worldness = zeros(1,total_num);

for p_temp = 1:total_num

    arr_small_worldness(p_temp) = gm_cell{p_temp}.small_worldness;

end

global_measures.small_worldness.mean =  mean( arr_small_worldness(:));
global_measures.small_worldness.std =  std( arr_small_worldness(:));
global_measures.small_worldness.all =   arr_small_worldness(:) ;
