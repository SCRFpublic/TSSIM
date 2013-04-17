function [ev_grid]=read_earthvision_grid(filename_and_path,g);
% This function reads a grid output by Petrel in EarthVision format

% Written by: Holly Michael
% April 2008

% Load pathname to SGEMS output data file
% filename_and_path - pathname to data file
% g: grid

% 
% Return - grid as a matlab array

fid=fopen(filename_and_path,'r');
if(fid==-1)
    fprintf('Error opening file %s\n',filename_and_path);
    error('read_SGEMS_output: error opening file');
end

% % Read header
 l=1;
 while l<=20       
     tline = fgetl(fid);    
     l=l+1;    
 end

% Read data
data=fscanf(fid,'%g',[1,inf]);
size(data)


fclose(fid);

x(:,1)=data(1,:);

for i=1:g.nx
    for j=1:g.ny
        %sim_grid(j,i)=x(i+g.nx*(j-1),1);
        ev_grid(g.ny-j+1,i)=x((i+g.nx*(j-1))*5-2,1);
        %ev_grid(j,i)=x((i+(j-1))*5-2,1);
    end
end
size(ev_grid)

figure;
%clims=[0 .1e-3];
imagesc(g.x,g.y,ev_grid);
title('Initial Topographic Elevation [m]');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
axis xy;
axis equal
