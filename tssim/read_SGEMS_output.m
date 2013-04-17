function [sim_grid]=read_SGEMS_output(filename_and_path,g,verbose);
% Read in SGEMS simulation
% Load pathname to SGEMS output data file
% fname - pathname to data file
% verbose - 1: comments otherwise nothing
% 
% Return - last SGEMS simulation
% x -data
% January 2008

% The original version of this file worked just as well.

% switch it around so it can be read by matlab (function written by Tapan Mukerji):
%[out, colnames, line1]=loadgeoeas(filename_and_path)

filename_and_path;
if ischar(filename_and_path)
fid3=fopen(filename_and_path);
line1=fgetl(fid3);
ncol=str2num(fgetl(fid3));
for k=1:ncol; colnames{k}=fgetl(fid3); end;
out=fscanf(fid3,'%f');
nrow=length(out)/ncol;
x=reshape(out,ncol,nrow).';
fclose(fid3);
end;


for i=1:g.nx
    for j=1:g.ny
        sim_grid(j,i)=x(i+g.nx*(j-1));
    end
end

if(verbose==1)
    fprintf('read_SGEMS_output: Read in last SGEMS simulation\n',filename_and_path,length(x));
end


