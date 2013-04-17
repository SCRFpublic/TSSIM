%retrieves information regarding core in SA.
% Inputs:
%     map: Conditional simulation area.
%     hard_data_location: matrix with hard data from get_harddata.
%     dirname: code locations obtained using get_dirname.
%     g: grid dimensions.
% Outputs:
%     length_core: length of core.
%     core_type: type: sand or shale.
%     type_contact: currently not used.
%     well_last_data: number of data intervals. If it's zero, all the well 
% conditions have been met.
%     datai: i-coord. of the well.
%     dataj: j-coord. of the well.
%     pos_data: position on the well. Sub2ind.

% written by Alejandro D. Leiva, June '09.
function [length_core core_type type_contact well_last_data datai dataj pos_data] = get_core_data_in_SA(map, hard_data_location,dirname,g)
pos_data = find(map.*hard_data_location(:,:,1) ~= 0);
pos_data = pos_data(1);
[datai dataj] = ind2sub([g.ny g.nx],pos_data);
well_number = hard_data_location(g.ny*g.nx + pos_data);
well_current_data = hard_data_location(2*g.ny*g.nx + pos_data);
well_last_data = hard_data_location(3*g.ny*g.nx + pos_data);
%reading the data
well_name = strcat(dirname,'\Well_data\w',int2str(well_number),'.txt');
well = importdata(well_name);
length_core = well(well_current_data+2,1);
core_type = well(well_current_data+2,2);
type_contact = well(well_current_data+1,3);