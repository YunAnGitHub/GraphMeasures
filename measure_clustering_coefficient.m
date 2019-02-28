function [ Clustering_Coefficient ] = measure_clustering_coefficient( gPPI_weighted_matrix)

% 2018-06-28 Yun-An Huang
% the script is used to claculate the clustering coefficent by using the
% formula introduced by Miyajima and sakuragawa 2014
% gPPI_weighted_matrix is a directed and weighted matrix.

% parameter setup

N_node = size(gPPI_weighted_matrix,1);
D_matrix = 1./gPPI_weighted_matrix;
Min_d = 1/max(gPPI_weighted_matrix(:));
C_array = zeros(1,N_node);

for ktemp = 1:N_node
    
    numerator = 0;
    denominator = 0;
    
    for itemp = 1:N_node
       
        for jtemp = 1:N_node

            numerator = numerator + 2/(0.5/(D_matrix(ktemp,itemp)+D_matrix(ktemp,jtemp))+D_matrix(itemp,jtemp));
            denominator = denominator + 2/(0.5/(D_matrix(ktemp,itemp)+D_matrix(ktemp,jtemp))+Min_d);
            
        end
    end
    
    
    C_array(ktemp) = numerator/denominator;
end

Clustering_Coefficient = sum(C_array)/N_node;

end