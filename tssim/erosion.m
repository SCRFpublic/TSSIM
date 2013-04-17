% Function that computes the erosion considering the gradient magnitud, 
% the gradient direction and the curvature profile.
% Inputs:
%     sim_layer: 2D geobody.
%     profc: Gradient direction map. Values 0-8.
%     grad: gradient magnitude map.
%     cond_thickness_map: lobe thickness map.
%     min_dir:
%     center_point: sediment source.
%     g: grid dimensions.
%     b: values used in the trend map constrction (see trend_map_maker)
% Output:
%     er: 2D erosion map ready to be stacked on base topo.
% Note: this function should be modified in case of conditional simulation 
% accounting for erosion.
% Wrritten by Alejandro D. Leiva
function er = erosion(sim_layer,profc,Grad,cond_thickness_map,min_dir,center_point,g,b)
eros_locs = find(sim_layer>0);
trend_map = trend_map_maker(g,center_point,b);

Grad = flipud(Grad);
min_g = min(min(Grad(eros_locs)));
Grad = Grad - min_g;
max_g = max(max(Grad(eros_locs)));
Grad_norm = Grad/max_g;
Grad_norm(sim_layer==0) = 0;
%Grad_norm(find(sim_layer==0)) = Grad_norm(find(sim_layer==0))*0.6;

cur_weigh = (profc>0)*0.8 + (profc<0)*1 + (profc==0)*0.9;
%cur_weigh(find(sim_layer==0)) = cur_weigh(find(sim_layer==0))*0.6;

dir_power = (min_dir==0)*1 + (min_dir==1)*0.8 +(min_dir==2)*0.6;
%dir_power(find(sim_layer==0)) = dir_power(find(sim_layer==0))*0.6;

ero = flipud(cur_weigh.*dir_power).*trend_map.*Grad_norm;
er_thick = ero/max(max(ero));

%erosion component due to mass weight
er = 0.2*cond_thickness_map/max(max(cond_thickness_map))+0.8*er_thick;
er = er/max(max(er));


