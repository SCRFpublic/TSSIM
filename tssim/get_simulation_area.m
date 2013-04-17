%function that gives simulations area.
% Inputs:
%     E: topography previously smoothened in case it is required.
%     pfield_startpoint: anchor point p-field.
%     show_plots: flag that indicates if map has to be plotted.
%     g: structure with grid dimensions.
%     min_river_length: condition on influence area length.
% Outputs:
%     sim_area: conditional simulation area.
%     start_point: anchor point drawn.
%     I: influence map for start_point.
%     pos_i: i-coord. anchor point.
%     pos_j: j-coord. anchor point.
% Written by Alejandro D. Leiva, June '09.
function [sim_area start_point I pos_i pos_j] = get_simulation_area(E,pfield_startpoint,show_plots,g,min_river_length)  
[M, N] = size(E);
R = dem_flow(E);
T = flow_matrix(E, R);
while (true) 
    [ri, rj]=get_new_start_point_indep(pfield_startpoint,g);
    rhs = zeros(numel(E), 1);
    idx = (rj-1)*M + ri;
    rhs(idx) = 1;
    I = T \ rhs;
    I = reshape(I, M, N);
    [riv_l pos_i pos_j] = river_length(I,[rj ri]);
    if riv_l>min_river_length && pos_i>10
        disp('Looking for Simulation Area....');
        A=zeros(M,N);
        A(I>0.01)=1;
        D = dependence_map(E, T, A);
        sim_area=zeros(M,N);
        sim_area(D>(prctile(D(D>0),25))) = 1;
        if show_plots == true
            AA=logical(A);
            vis_map(D, E,AA);
            disp('done!');
            title('Influence and Dependence Map')
        end
        break;
    end
end
start_point=[rj ri];