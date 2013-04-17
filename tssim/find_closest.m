%area is a image with zeros (non-area) and ones (area of interest)
%the area might be a perimeter as well
%point = [i,j]
function [I J] = find_closest(area,point)
    [M N] = size(area);
    index = find(area==1);
    [INDI INDJ] = ind2sub([M N],index);
    dist = sqrt((INDI-point(2)).^2+(INDJ-point(1)).^2);
    min_position = find(dist == min(dist));
    I = INDI(min_position);
    J = INDJ(min_position);
  