function [ dist_array, track_cell ] = dijkstra_weighted( dist_matrix, source )
%DIJKSTRA_WEIGHTED Summary of this function goes here
%   Detailed explanation goes here

% 2018-Mar-20 Yun-An Huang
% track for the shortest path with multi pathway.

% 2018-feb-06 Yun-An Huang
% the implementation of dijkstra algorithm

% input:
% 1. dist_matrix: the n x n matrix store the distance of edge, raw is
% output, column is input.
% 2. source: the source

% output:
% 1. dist_array: the distant array and the track with cell structure.
% 2. track_cell: the track of different nodes connected to.


%% initialization
node_num = size(dist_matrix,1); % the node number

Unreach_node = ones(1,node_num); % the node haven't been reached.

dist_array = inf(1,node_num); % initialize the dist_array. 
dist_array(source) = 0; % the distant of self node is zero.

% initialize the array of the previous node. by tracking the previous node we can construct the tracks from source to other destination.
for itemp = 1:node_num
    prev_array{itemp} = 0; 
end

%% start to search the minimum path length of each node.

while sum(Unreach_node)~=0
    current_node_temp = find(dist_array==min(dist_array(Unreach_node>0)) & Unreach_node>0); % select the least distant node as the current node.
    current_node= current_node_temp(1);
    
    Unreach_node(current_node) = 0; % the node is reached.
    neighbor_of_current_node = find(dist_matrix(current_node, : )<inf & Unreach_node >0); % find the neighbor of current node
    
    for itemp = neighbor_of_current_node
        % for every neighbor update the distant
        new_dist = dist_array(current_node) + dist_matrix(current_node,itemp);
        if new_dist < dist_array(itemp)
           
            dist_array(itemp) = new_dist; % update the path distant.
            prev_array{itemp} = current_node; % update the previous node
            
        elseif new_dist == dist_array(itemp)
            
            prev_array{itemp} = [prev_array{itemp} current_node]; % the multiple track
        end
        
    end
            
end


%% track the path for each node.

track_cell =[];

for itemp = 1:node_num % track for each node.
    

    track_num = 1; % the number of tracks.
    track_temp = 1;
   
    % initialize
    current_node = itemp;
    track_array = [current_node];
    track_cell{itemp,track_temp}=track_array;
    
    
    while track_temp<=track_num % track for multiple tracks

        track_array = track_cell{itemp,track_temp};
        current_node = track_array(1);
        
        while prev_array{current_node}~=0 % iteratively track the node
    
                
            prev_array_temp = prev_array{current_node};
            for p_temp = 1:length(prev_array_temp) 

                if p_temp ==1
                    track_array_temp = [prev_array_temp(p_temp) track_array]; % add the prev node to the track
                    track_cell{itemp,track_temp} = track_array_temp;

                else
                    track_array_temp = [prev_array_temp(p_temp) track_array];
                    track_cell{itemp,track_num+p_temp-1} = track_array_temp; % the new track
                end

            end

            track_num = track_num+length(prev_array_temp)-1;
            
            track_array=track_cell{itemp,track_temp};
            current_node = track_array(1);
            
        end
    
        track_temp = track_temp +1;
   
    end
    
    
end




end

