function [max_dist pos_i pos_j]= river_length(I,start_point)
%Used by tssim to compute river length

[i j] = ind2sub(size(I),find(I>0.95));
dist = sqrt((i-start_point(2)).^2+(j-start_point(1)).^2);
max_dist = max(dist);
pos_i = i(dist==max_dist);
pos_j = j(dist==max_dist);
pos_i = pos_i(1);
pos_j = pos_j(1);