% write a grid in GEOEAS/GSLIB format to a named file
% write_gslib_grid (fname,data,g,dim)
% fname - file name
% data - data matrix
% g - structure containing grid dimension data
% verbose - 1: comments otherwise nothing

% Written by: Jonathan Steward/Holly Michael

function write_gslib_grid (fname,data,g,dim)
tic;

fid=fopen(fname,'w');
if(fid==-1),
    fprintf('write_xydat_grid: Error opening file %s\n',fname);
    error('cannot open file');
end

fprintf(fid,'Description: Lobe Thickness Data\n');

if dim==3
    fprintf(fid,'%d\n',4);
    fprintf(fid,'Field 1: x\n');
    fprintf(fid,'Field 2: y\n');
    fprintf(fid,'Field 3: layer\n');
    fprintf(fid,'Field 4: thickness or lobe\n');
    
    for k=1:size(data,3)
         for j=1:g.ny
             for i=1:g.nx
                 fprintf(fid,'%.6f %.6f %.6f %.6f\n',g.x(i),g.y(j),k,data(j,i,size(data,3)-k+1));
             end
         end
    end

elseif dim==2
    fprintf(fid,'%d\n',3);
    fprintf(fid,'Field_1:_x\n');
    fprintf(fid,'Field_2:_y\n');
    fprintf(fid,'Field_3:_thickness_or_lobe\n');
    
     for j=1:g.ny
        for i=1:g.nx
         fprintf(fid,'%.6f %.6f %.6f\n',g.x(i),g.y(j),data(j,i));
        end
     end
else
       error('not a valid dimension');  
end


fclose(fid);
