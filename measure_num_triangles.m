function [ num_triangles_arr ] = measure_num_triangles( gPPI_matrix )
%MEASURE_NUM_TRIANGLES Summary of this function goes here
%   Detailed explanation goes here


% Yun-An Huang 2018-March-26
% calculate the number of triangle of each nodes of a gPPI matrix.
% T(i) = 1/2*sum( (w(i,j)^(1/3)+w(j,i)^(1/3))*(w(i,h)^(1/3)+w(h,i)^(1/3))*(w(j,h)^(1/3)+w(h,j)^(1/3)))) for
% j,h is the arbitary node.

% input:
% 1. gPPI_matrix: store the weighted of gPPI matrix.

% output:
% 1. num_triangles_arr: the array store the number of triangles of each
% nodes.



node_num = size(gPPI_matrix,1);
num_triangles_arr = zeros(1,node_num);

% calculate triangles of each nodes.
for itemp = 1:node_num

    triangle_temp_matrix = zeros(node_num,node_num);
    
    node_temp_array = 1:node_num;
    node_temp_array(itemp)=[]; % do not consider the source node.
    
    % calculate the triangles between i and arbitary j,k
    for jtemp = node_temp_array
        for ktemp = node_temp_array
       
            if jtemp == ktemp
            
                triangle_temp_matrix(jtemp,ktemp)=0;
                
            else
                triangle_temp_matrix(jtemp,ktemp)=((gPPI_matrix(itemp,jtemp)^(1/3)+gPPI_matrix(jtemp,itemp)^(1/3))*(gPPI_matrix(itemp,ktemp)^(1/3)+gPPI_matrix(ktemp,itemp)^(1/3))*(gPPI_matrix(jtemp,ktemp)^(1/3)+gPPI_matrix(ktemp,jtemp)^(1/3)));
                
            end
                
        
        end
    end
    
   
    num_triangles_arr(itemp)=sum(triangle_temp_matrix(:))/2;
    
end



end

