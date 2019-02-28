function [ betweenness_centrality ] = measure_betweeness_centrality( track_shortpath, node_i )
%MEASURE_BETWEENESS_CENTRALITY Summary of this function goes here
%   Detailed explanation goes here

% 2018-March-21 Yun-An Huang
% this function is used to calculate the betweenness centrality.

% input:
% 1. the track cell contained track paths from global measures. (track_shortpath)
% 2. the identifity node.

% output:
% 1. the betweenness centrality of node i.


%% initialize
node_num = size(track_shortpath,1); % the matrix size

betweenness_count = 0;

%% track for each path

for itemp = 1:node_num
    for jtemp = 1:node_num
 
        if itemp~= node_i && jtemp ~= node_i && itemp ~=jtemp 
        
            count_all = 0; % all tracks
            count_i = 0; % the tracks pass through node i
            
            track_cell = track_shortpath{itemp,jtemp};
            for t_temp  = 1:length(track_cell)
                
                if isempty(track_cell{t_temp}) % nothing in the track
        
                else
                    count_all = count_all+1;
                    if sum(track_cell{t_temp}==node_i) % the track pass through node_o
                        count_i = count_i + 1;
                     
                    end
                    
                end
                
            end
            
            
            betweenness_count = betweenness_count+count_i/count_all ;
            
        end
        
        
        
    end
end

betweenness_centrality = betweenness_count/(node_num-1)/(node_num-2);








end

