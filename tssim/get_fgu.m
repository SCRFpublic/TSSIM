% Gets fine-grained unit considering an altitude component. Thickner at low 
% locations.
% inputs:
%     prev_sim_elev: Base topography.
%     fines_thick: Average thickness drawn from CDF.
%     altitude_component: [0,1] value that tells how much the altitude of 
% the previous topography has to be considered in the sim.
% Output:
%     fgu: 2D matrix that is to be stacked as IFGU.
% By Alejandro D. Leiva, June '09.

function fgu = get_fgu(prev_sim_elev,fines_thick,altitude_component)
[m n] = size(prev_sim_elev);
fgu = ones(m,n)*(1-altitude_component)*fines_thick;
max_elev = max(max(prev_sim_elev));
weight = abs(max_elev - prev_sim_elev)/max_elev;
fgu = fgu + weight*fines_thick*altitude_component;