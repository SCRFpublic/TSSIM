% Writes .par file for sgsim to add noise to the simulation. 
% Only for non-conditional simulation. Erosion
% Inputs:
%     filename: file name desired. 
%     gx: x-grid size.
%     gy: y-grid size.
% Outputs:
%     filename: .par file 

% Written by Alejandro D. Leiva, June '09.
function write_dimension_parfile_erosion(filename,gx,gy)
fid =  fopen(filename,'w');
fprintf(fid,'                  Parameters for SGSIM\r\n');
fprintf(fid,'                  ********************\r\n');
fprintf(fid,'                                        \r\n');
fprintf(fid,'START OF PARAMETERS:\r\n');
fprintf(fid,'hard_data_erosion.txt                            -file with data\r\n');
fprintf(fid,'1  2  3  4  0  0              -  columns for X,Y,Z,vr,wt,sec.var.\r\n');
fprintf(fid,'-1.0       1.0e21             -  trimming limits\r\n');
fprintf(fid,'1                             -transform the data (0=no, 1=yes)\r\n');
fprintf(fid,'sgsim.trn                     -  file for output trans table\r\n');
fprintf(fid,'0                             -  consider ref. dist (0=no, 1=yes)\r\n');
fprintf(fid,'hist.txt                  -  file with ref. dist distribution\r\n');
fprintf(fid,'1  0                          -  columns for vr and wt\r\n');
fprintf(fid,'0    1.0                   -  zmin,zmax(tail extrapolation)\r\n');
fprintf(fid,'1       1.0                   -  lower tail option, parameter\r\n');
fprintf(fid,'1      2.0                   -  upper tail option, parameter\r\n');
fprintf(fid,'1                             -debugging level: 0,1,2,3\r\n');
fprintf(fid,'sgsim.dbg                     -file for debugging output\r\n');
fprintf(fid,'sgs_erosion.txt                     -file for simulation output\r\n');
fprintf(fid,'1                             -number of realizations to generate\r\n');
fprintf(fid,strcat(num2str(gx),'    1.0    1.0              -nx,xmn,xsiz\r\n'));
fprintf(fid,strcat(num2str(gy),'    1.0    1.0              -ny,ymn,ysiz\r\n'));
fprintf(fid,'1     0.5    1.0              -nz,zmn,zsiz\r\n');
fprintf(fid,strcat(num2str(ceil(rand(1)*10000)),'                         -random number seed\r\n'));
fprintf(fid,'2     12                       -min and max original data for sim\r\n');
fprintf(fid,'10                            -number of simulated nodes to use\r\n');
fprintf(fid,'1                             -assign data to nodes (0=no, 1=yes)\r\n');
fprintf(fid,'1     3                       -multiple grid search (0=no, 1=yes),num\r\n');
fprintf(fid,'0                             -maximum data per octant (0=not used)\r\n');
fprintf(fid,'60.0  60.0  60.0              -maximum search radii (hmax,hmin,vert)\r\n');
fprintf(fid,'0.0   0.0   0.0              -angles for search ellipsoid\r\n');
fprintf(fid,'51    51    11                -size of covariance lookup table\r\n');
fprintf(fid,'1     0.0   0.0              -ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC\r\n');
fprintf(fid,'local_varying_mean.txt             -  file with LVM, EXDR, or COLC variable\r\n');
fprintf(fid,'1                             -  column for secondary variable\r\n');
fprintf(fid,'1    0.0                      -nst, nugget effect\r\n');
fprintf(fid,'1    1  0.0   0.0   0.0     -it,cc,ang1,ang2,ang3\r\n');
fprintf(fid,'         100.0  100.0  1.0     -a_hmax, a_hmin, a_vert1\r\n');
fclose(fid)