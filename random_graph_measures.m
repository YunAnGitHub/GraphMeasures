function [ global_measures ] = random_graph_measures( gPPI_weighted_matrix, permute_num )
%GPPI_GRAPH_MEASURES Summary of this function goes here
%   Detailed explanation goes here

% Yun-An Huang 2018-Apr-17
% to randomize the non-diagonal entries.

% Yun-An Huang 2018-Mar-02
% this function is calculating the global measures with the randomize
% network with same data distribution.
%

rng('shuffle');

gm={};

node_num = size(gPPI_weighted_matrix,1);
parfor p_temp = 1:permute_num

    
    
    
    
%randomize with diagonal entris
%     rand_ix = reshape(randperm(length(gPPI_weighted_matrix(:))),size(gPPI_weighted_matrix));
%     gPPI_weighted_matrix_permute = gPPI_weighted_matrix(rand_ix);

%randomize with non-diagonal entris

    idx = eye(node_num);
    gPPI_weighted_matrix_permute = zeros(node_num,node_num);
    
    non_diag_data  = gPPI_weighted_matrix(~idx);
    idx_rd = randperm(length(non_diag_data));
    gPPI_weighted_matrix_permute(~idx) = non_diag_data(idx_rd);
    
    gm{p_temp}=gPPI_weighted_graph_measures(gPPI_weighted_matrix_permute);
           
end



global_measures = mean_gm(gm);

