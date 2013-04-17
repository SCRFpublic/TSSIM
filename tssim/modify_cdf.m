function modify_cdf(file_name,end_point)
% Used by tssim to modify cdf by extrapolating endpoints.
    fid =  fopen(file_name,'w');
    fprintf(fid,'%6.2f %6.2f\r\n',end_point,1);
    fprintf(fid,'%6.2f %6.2f\r\n',0,0);
    fclose(fid)
end