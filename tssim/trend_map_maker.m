%Program generates trend map from sediment source until b.
%Input:
%       g : structure with grid dimensions.
%       b<=1 : values are zeros at b of the sortest side of the X-Y model.
%Output:
%       trend_map_normavg: trend map starting from sediment source until b.
%Written by H. Michael.
%Modified by A. Leiva, June '09.
function trend_map_normavg=trend_map_maker(g,center_point,b)
sed_source_loc = center_point;
trend_map=ones(g.ny,g.nx,1);
max_dist_cells=round(b*min(g.nx,g.ny));   
for iii=1:g.nx
    for jjj=1:g.ny
        dx_cells=abs(iii-sed_source_loc(1));
        dy_cells=abs(jjj-sed_source_loc(2));
        dist_source=sqrt(dx_cells^2+dy_cells^2);
        if dist_source/max_dist_cells<1 
            trend_map(jjj,iii,1)=1-dist_source/max_dist_cells;
        else
            trend_map(jjj,iii,1)=0;
        end
    end
end
avg_trend_map=mean(mean(trend_map));
trend_map_normavg=trend_map/avg_trend_map;