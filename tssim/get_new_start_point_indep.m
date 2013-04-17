function [row_pick, col_pick]=get_new_start_point_indep(pmap_combined_norm,g)

%This function uses a probability map to get the next lobe start point
% In this version, each new point is independent of the
% previously-simulated point.

% Written by: Holly Michael
% April 2008


% Use a Cox's point process to randomly pick a point based on the probability map:
    % Randomly choose a cell. Then retain the cell based on its
    % probability. If not retained, randomly choose another cell and repeat
    % (check that this is really a Cox Process (perhaps similar to a
    % poisson process?).
    keep=0;
    count=0;
    while keep==0
        % Random Pick for i,j:
        ipick=ceil(g.nx*rand);
        jpick=ceil(g.ny*rand);
      %  pmap_combined_norm(jpick,ipick);
        keep=pmap_combined_norm(jpick,ipick)>=rand;
        count=count+1;
    end
    row_pick=jpick;
    col_pick=ipick;
    count;
   