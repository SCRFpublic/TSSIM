function g=make_grid(g,verbose)
% Calculate grid dimensions
% g - structure containing grid origin, length and spacing
% verbose - 1: comments otherwise nothing
% Return
% g - structure containing calculated grid information

%Written by Jonathan Stewart
%Edited by Hongmei Li ---extend 2D grid to 3D

g.xmax=g.x0+g.lx;
g.x=g.x0+g.dx/2:g.dx:g.xmax-g.dx/2;
g.nx=length(g.x);

g.ymax=g.y0+g.ly;
g.y=g.y0+g.dy/2:g.dy:g.ymax-g.dy/2;
g.ny=length(g.y);

%g.zmax=g.z0+g.lz;
%g.z=g.z0+g.dz/2:g.dz:g.zmax-g.dz/2;
%g.nz=length(g.z);

g.nxyz=g.nx*g.ny; %*g.nz;

if(verbose)
   fprintf('make_grid: nx %d ny %d\n',g.nx,g.ny);
   fprintf('make_grid: xmin %.2f xmax %.2f ymin %.2f ymax %.2f\n',g.x0,g.xmax,g.y0,g.ymax);
   fprintf('make_grid: dx %.2f dy %.2f\n',g.dx,g.dy);   
   fprintf('make_grid: x(1) %.2f x(nx) %.2f\n',g.x(1),g.x(g.nx));
   fprintf('make_grid: y(1) %.2f y(ny) %.2f\n',g.y(1),g.y(g.ny));
   %fprintf('make_grid: z(1) %.2f z(nz) %.2f\n',g.z(1),g.z(g.nz));
end
