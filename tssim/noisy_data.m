% Writes noisy data called for .par file
% Inputs:
% filename
% n_points_fraction ( 0 <=x<=1), ref.: 0.5
% by Alejandro D. Leiva, june '09
function noisy_data(filename,thickness_map,n_points_fraction,g)
random_visit=randperm(g.nx*g.ny);
fid =  fopen(filename,'w');
fprintf(fid,'data\r\n');
fprintf(fid,'%.0f\r\n',4);
fprintf(fid,'x\r\n');
fprintf(fid,'y\r\n');
fprintf(fid,'z\r\n');
fprintf(fid,'data\r\n');
for i = 1:75000*n_points_fraction
    [y,x] = ind2sub([g.ny g.nx],random_visit(i));
    w = [x y 0.5 thickness_map(random_visit(i))];
    fprintf(fid,'%6.2f  %6.2f %6.2f  %6.2f\r\n',w);
end
fclose(fid)