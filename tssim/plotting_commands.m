% function plotting_commands
% 
% 
% 
% 
% eros_locs=find(sim_layer>0);
%        [profc,planc]=curvature(prev_sim_elev);
%         surf(prev_sim_elev,profc); title('Profile curvature')
%         
%         
%         
%         
% surf(flipud(prev_sim_elev,profc));set(gca,'DataAspectRatio',[1 1 0.3],'PlotBoxAspectRatio',[1 1 0.3])
% shading interp; camlight
% 
% [DEM] = base_topo;
% [G,ASP] = gradient8(DEM,1);
% surf(DEM,G); shading interp; camlight

h = surf(surface_cube(:,:,1));set(gca,'DataAspectRatio',[1 1 0.3],'PlotBoxAspectRatio',[1 1 0.3]);
shading interp; camlight;
rotate(h,[0 0 1],180);
