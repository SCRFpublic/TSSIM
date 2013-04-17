% gets angle used in the object-based simulation in the lobe orientation. The 
% direction is in the highest prc-percetile of furthest point on the SA
% boundary to the anchor point (start_point).
% Inputs:
%     start_point: anchor point.
%     countour_map: simulation area contour map.
%     prc: percentile used in the point selection on the boundary.
% Ouput:
%     angle: measured from east direction clock-wise [0-180].
%by Alejandro Leiva, June '09
function angle = get_angle(start_point,contour_map,prc)
%start_point(i,j)
index = find(contour_map==1);
[m n] = size(contour_map);
[IND1 IND2] = ind2sub([m n],index);
dist = sqrt((IND1-start_point(1)).^2+(IND2-start_point(2)).^2);
large_values = find(dist>prctile(dist,prc));
dy = -start_point(1)+IND1(large_values);
dx = -start_point(2)+IND2(large_values);
angle = 180*atan2(dy,dx)/pi;
end