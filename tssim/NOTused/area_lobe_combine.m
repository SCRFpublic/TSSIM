function [thickness_map closeBW] = area_lobe_combine(sim_area_bi,lobe_ind_map,multiplier,start_point,chan_ind,g)
%used internally by tssim to combine simulation area and lobe area.

w = sim_area_bi+3*lobe_ind_map;
ww=zeros(g.ny,g.nx,2);
ww(:,:,1)=w;
ww_smooth=smooth3(ww,'box',[31 31 1]);
ww_sm=ww_smooth(:,:,1);
lobe_ind_map=ww_sm>ww_sm(start_point(2),start_point(1))*multiplier;
ind_map=lobe_ind_map+fliplr(flipud(chan_ind));
se = strel('disk',7);
closeBW = imclose(ind_map>0,se);
%% Lobe Thickeness Assignment
thickness_map = lobe_thickness_assignment(closeBW,start_point,g);
%figure;imagesc(thickness_map);colorbar
