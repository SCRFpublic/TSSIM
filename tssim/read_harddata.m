% reads hard data from text files in tssim\Well_data\
% Inputs:
%     conditional: conditional flag.
%     n_wells: number of wells.
%     dirname: folder where the code is located. Easily obtained using 
% get_dirname.
%     g: grid specifications.
% Outputs:
%     conditional_interval: number of data conditioning the simulation.
%     hard data location: matrix of x-size. First X-Y map has well number, 
%second the current position and the third the last data. 
%     warping_flag: flag that tells if grid deformation is necessary.

% by Alejandro Leiva
    
function [conditional_interval hard_data_location warping_flag]= ...
    read_harddata(conditional,n_wells,dirname,g)
hard_data_location = zeros(g.ny,g.nx,4);
if conditional == true
    warping_flag = false;
    for i=1:n_wells
        well_name = strcat(dirname,'\Well_data\w',int2str(i),'.txt');
        well = importdata(well_name);
        hard_data_location(well(1,2),well(1,1),1) = 1;
        hard_data_location(well(1,2),well(1,1),2) = i;%well number
        hard_data_location(well(1,2),well(1,1),3) = 0;%current data position
        hard_data_location(well(1,2),well(1,1),4) = size(well,1)-1;%last data
    end
    conditional_interval = sum(sum(hard_data_location(:,:,4)));  
else
    conditional_interval = 0;
    warping_flag = false;
end