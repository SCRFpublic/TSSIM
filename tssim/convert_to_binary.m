function convert_to_binary(fname,input_path_name, output_path_name,g2s_code_path, ML_code_path,dtype,gridname, g)

% This function converts an ascii TI file to binary to be used in SGEMS
%   Written by Holly Michael
%   January, 2008

% Inputs:
%  fname: the filename of the g2s parameter file
%  input_path_name: full path and filename for gslib format file
%  output_path_name: full path and filename for binary file (will be created)
%  g2s_code_path: path to exe file for g2s
%  matlab_code_path: path to matlab code
%  type: 0 for pointset, 1 for grid
%  gridname: name of the grid

cd(g2s_code_path)
% Write the parameter file for geoeas2sgems:
fid=fopen(fname,'wt');

if(fid==-1),
    fprintf('write_xydat_grid: Error opening file %s\n',fname);
    error('cannot open file');
end

fprintf(fid,'%s\n',input_path_name);
fprintf(fid,'%s\n',output_path_name);
fprintf(fid, '%s\n',gridname);
fprintf(fid, '%d\t %s\n',dtype,'%% type of object: pointset (0), grid (1)');
fprintf(fid, '%s\t %s\n','0', '%% If pointset: 2D (0) or 3D (1)');
fprintf(fid, '%d %d %d\n', g.nx,g.ny,1);
fprintf(fid, '%d %d %d\n', 1,1,1);
fprintf(fid, '%d %d %d\n', 0,0,0);
fprintf(fid, '%d\n', 1);
fprintf(fid, '%d\n', -99999);

fclose(fid);

system('geoeas2sgems')
cd(ML_code_path)
