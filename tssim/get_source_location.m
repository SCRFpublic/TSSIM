%Finds sediment source location in the range.
%Imputs: 
%   prev_sim_elev: current simulations topography.
%   range: range value for the sediment source location on the upper
%          boundary of the model
%   start_point: anchor_point.
%Output:
%   s_source_location = sediment source location (x,y)
%by Alejandro Leiva, June '09
function s_source_location = get_source_location(prev_sim_elev,start_point,range)
win = ceil(range/2);
[m n] = size(prev_sim_elev);
if start_point <= floor(n/2)
    seg = prev_sim_elev(1,floor(n/2)-win:floor(n/2));
    pos = find(seg == min(seg));
    s_source_location = [floor(n/2)-win+pos(1) 1];
else
    seg = prev_sim_elev(1,floor(n/2)+1:floor(n/2)+win);
    pos = find(seg == min(seg));
    s_source_location = [floor(n/2)+1+pos(1) 1];
end
    
    
    