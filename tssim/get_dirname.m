% Returns fullpath

% Written by Alejandro D. Leiva
% June, 2009
function dirname = get_dirname;
fp = mfilename('fullpath');
dirname = fileparts(fp);
slash = strfind(dirname,'\');
dirname = dirname(1:slash(end)-1);