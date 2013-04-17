%Returns river flow direction according to the following pattern:
%                    8 1 2
%                    7 c 3
%                    6 5 4
% Inputs:
%     I: Influence map.
%     start_point: anchor point (x,y).
% Outputs:     
%     r_dir: river direction (1-8)

% Written by Alejandro D. Leiva.

function r_dir = river_direction(I,start_point)
[i j] = ind2sub(size(I),find(I>0.95));
dist = sqrt((i-start_point(2)).^2+(j-start_point(1)).^2);
max_dist_index = find(dist==max(dist));
angle = atan((i(max_dist_index)-start_point(2))/(j(max_dist_index)-start_point(1)));
if -pi/2<=angle && angle<=-3*pi/8
    r_dir = [1 5];
elseif -3*pi/8<=angle && angle<=-pi/8
    r_dir = [2 6];
elseif -pi/8<=angle && angle<=pi/8
    r_dir = [3 7];
elseif pi/8<=angle && angle<=3*pi/8
    r_dir = [4 8];
elseif 3*pi/8<=angle && angle<=pi/2
    r_dir = [5 1];
end
    