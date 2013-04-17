function [x]=read_cdf_data(fname,verbose);
% Load pattern parameter cdf data file
% fname - name of file
% verbose - 1: comments otherwise nothing
% 
% Return - 2 columns of cdf file
% x -data
% y- percentage
% n - number of lines

% Written by: Hongmei Li
% Date: 2007 Summer

fid=fopen(fname,'r');
if(fid==-1)
    fprintf('Error opening file %s\n',fname);
    error('read_cdf: error opening file');
end

% Read header
% l=1;
% while l<=5       
%     tline = fgetl(fid);    
%     l=l+1;    
% end

% Read data
data=fscanf(fid,'%g %g',[2,inf]);

fclose(fid);

x(:,1)=data(1,:);
x(:,2)=data(2,:);

if(verbose==1)
    fprintf('read_cdf: %s: %d points\n',fname,length(x));
end

