function [ global_measures ] = gPPI_weighted_graph_measures_normalize( gPPI_weighted_matrix , permute_num)
%GPPI_GRAPH_MEASURES Summary of this function goes here
%   Detailed explanation goes here

% 2018-Apr-18 Yun-An Huang
% calculate the percentage of permutation.

% Yun-An Huang 2018-March-02
% this function is calculating the normalized global measures.
% the global meausures is normalize by the random network with the same
% data distribution.


gm_ori = gPPI_weighted_graph_measures(gPPI_weighted_matrix);
gm_rand = random_graph_measures(gPPI_weighted_matrix,[permute_num]);


%% Sum of weights


global_measures.sum_of_weight = gm_ori.sum_of_weight./gm_rand.sum_of_weight.mean;
global_measures.all.sum_of_weight = gm_rand.sum_of_weight.all;
global_measures.ori.sum_of_weight = gm_ori.sum_of_weight;
% global_measures.percentage.sum_of_weight = sum(gm_rand.sum_of_weight.all>= gm_ori.sum_of_weight)/length(gm_rand.sum_of_weight.all);
% the percentage of sum_of_weight is meaningless, since the distribution is
% the same.

%% characteristic path length.

global_measures.char_path_length = gm_ori.char_path_length./gm_rand.char_path_length.mean;
global_measures.percentage.char_path_length = sum(gm_rand.char_path_length.all>= gm_ori.char_path_length)/length(gm_rand.char_path_length.all);
global_measures.all.char_path_length = gm_rand.char_path_length.all;
global_measures.ori.char_path_length = gm_ori.char_path_length;

%% global efficiency

global_measures.global_eff = gm_ori.global_eff./gm_rand.global_eff.mean;
global_measures.percentage.global_eff = sum(gm_rand.global_eff.all>= gm_ori.global_eff)/length(gm_rand.global_eff.all);
global_measures.all.global_eff = gm_rand.global_eff.all;
global_measures.ori.global_eff = gm_ori.global_eff;

%% clustering coefficient.

global_measures.cluster_coeff = gm_ori.cluster_coeff./gm_rand.cluster_coeff.mean;
global_measures.percentage.cluster_coeff = sum(gm_rand.cluster_coeff.all>= gm_ori.cluster_coeff)/length(gm_rand.cluster_coeff.all);
global_measures.all.cluster_coeff = gm_rand.cluster_coeff.all;
global_measures.ori.cluster_coeff = gm_ori.cluster_coeff;

%% Transitivity

global_measures.transitivity = gm_ori.transitivity./gm_rand.transitivity.mean;
global_measures.percentage.transitivity = sum(gm_rand.transitivity.all>= gm_ori.transitivity)/length(gm_rand.transitivity.all);
global_measures.all.transitivity = gm_rand.transitivity.all;
global_measures.ori.transitivity = gm_ori.transitivity;

%% local efficiency

global_measures.local_efficiency = gm_ori.local_efficiency./gm_rand.local_efficiency.mean;
global_measures.percentage.local_efficiency  = sum(gm_rand.local_efficiency .all >= gm_ori.local_efficiency )/length(gm_rand.local_efficiency.all);
global_measures.all.local_efficiency  = gm_rand.local_efficiency.all;
global_measures.ori.local_efficiency = gm_ori.local_efficiency;

%% modularity
% since this is a weighted matrix without setting a threshold, every nodes
% is in the same module.


global_measures.modularity  = gm_ori.modularity./gm_rand.modularity.mean;
global_measures.percentage.modularity  = sum(gm_rand.modularity .all >= gm_ori.modularity )/length(gm_rand.modularity.all);
global_measures.all.modularity  = gm_rand.modularity.all;
global_measures.ori.modularity  = gm_ori.modularity;

%% assortativity coefficient


global_measures.assortativity  = gm_ori.assortativity./gm_rand.assortativity.mean;
global_measures.percentage.assortativity  = sum(gm_rand.assortativity.all >= gm_ori.assortativity )/length(gm_rand.assortativity.all);
global_measures.all.assortativity  = gm_rand.assortativity.all;
global_measures.ori.assortativity  = gm_ori.assortativity;

%% Small worldness

global_measures.small_worldness = gm_ori.small_worldness./gm_rand.small_worldness.mean;
global_measures.percentage.small_worldness  = sum(gm_rand.small_worldness.all >= gm_ori.small_worldness )/length(gm_rand.small_worldness.all);
global_measures.all.small_worldness  = gm_rand.small_worldness.all;
global_measures.ori.small_worldness = gm_ori.small_worldness;
