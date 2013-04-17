%function that gives simulations area that fits the hard data (data)
% Inputs:
%     E: topography previously smoothened in case it is required.
%     data: data coordinates.
%     pfield_startpoint: anchor point p-field.
%     g: structure with grid dimensions.
% Outputs:
%     sim_area: conditional simulation area.
%     start_point: anchor point drawn.
%     checked_points: points checked during iterations.
%     I: influence map for start_point.
%     pos_i: i-coord. anchor point.
%     pos_j: j-coord. anchor point.
% Written by Alejandro D. Leiva, June '09.    
function [sim_area start_point checked_points I pos_i pos_j] = get_simulation_area_conditional(E,data,pfield_startpoint,g)  
[M, N] = size(E);
R = dem_flow(E);
T = flow_matrix(E, R);
checked_points = zeros(M,N);
data_matrix = zeros(M,N);
data_matrix(data(:,2),data(:,1)) = 1;
while (true) 
    [ri, rj]=get_new_start_point_indep(pfield_startpoint,g);
    rhs = zeros(numel(E), 1);
    idx = (rj-1)*M + ri;
    rhs(idx) = 1;
    I = T \ rhs;
    I = reshape(I, M, N);
    [riv_l pos_i pos_j] = river_length(I,[rj ri]);
    if riv_l > 15 && checked_points(ri,rj)~=1
        disp('Looking for Conditional Simulation Area....');
        A=zeros(M,N);
        A(I>0.01)=1;
        D = dependence_map(E, T, A);
        sim_area=zeros(M,N);
        sim_area(D>0) = 1;
        checked_points(D>(prctile(D(D>0),10))) = 1;
        %figure;imagesc(checked_points);set(gcf,'Color',[1 1 1]);    
        if sum(sum(sim_area.*data_matrix~=0))
            disp('done');
            break;       
        end   
    end
end
start_point=[rj ri];