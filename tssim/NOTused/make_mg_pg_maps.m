function [prob_map_mg, prob_map_pg]=make_mg_pg_maps(mg_cdf, pg_cdf,sed_source_loc,g,verbose, plot_flag)
%MAKE_MG_PG_MAPS Outputs: 
% Outputs: 
%     prob_map_mg: map of probabilities at every point on the grid for migration based on the sediment source
%     prob_map_pg: map as above, for progradation

% Inputs:
%   mg_cdf: CDF of migration distances
%   pg_cdf: CDF of progradation distances
%   sed_source_loc: location in cells of sediment source
%   g: simulation grid
%   verbose: whether or not to print out progress (used in read_pattern_cdf
%   plot_flag: whether or not to plot the figures (1 for plot, 0 for no plot)

%read CDFs 
% Make a grid of progradation and migration probabilities
flip_mg_cdf=flipud(mg_cdf);
flip_pg_cdf=flipud(pg_cdf);
% Make a PDF from the CDF:
mg_pdf=zeros(length(mg_cdf(:,1))+1,2);
pg_pdf=zeros(length(pg_cdf(:,1))+1,2);
max_mgdist=max(mg_cdf(:,1));
max_pgdist=max(pg_cdf(:,1));

mg_pdf(1,1)=mg_cdf(end,1);  % first row is the distance, second is the prob.
pg_pdf(1,1)=pg_cdf(end,1);
mg_pdf(end,1)=max_mgdist;
pg_pdf(end,1)=max_pgdist;
for i=2:length(mg_cdf(:,1))
    mg_pdf(i,1)=flip_mg_cdf(i-1,1)+(flip_mg_cdf(i,1)-flip_mg_cdf(i-1,1))/2;
    mg_pdf(i,2)=(flip_mg_cdf(i,2)-flip_mg_cdf(i-1,2))/(flip_mg_cdf(i,1)-flip_mg_cdf(i-1,1));
end
for i=2:length(pg_cdf(:,1))
    pg_pdf(i,1)=flip_pg_cdf(i-1,1)+(flip_pg_cdf(i,1)-flip_pg_cdf(i-1,1))/2;
    pg_pdf(i,2)=(flip_pg_cdf(i,2)-flip_pg_cdf(i-1,2))/(flip_pg_cdf(i,1)-flip_pg_cdf(i-1,1));
end
min_mgdist=min(mg_pdf(:,1));
min_pgdist=min(pg_pdf(:,1));
for i=1:g.nx
    for j=1:g.ny
        dx=(i-sed_source_loc(1))*g.dx;
        dy=abs(j-sed_source_loc(2))*g.dy;
        if dx>=max_mgdist|dx<=min_mgdist
            prob_map_mg(j,i)=0;
        else
            prob_map_mg(j,i)=interp1(mg_pdf(:,1), mg_pdf(:,2),dx);
        end
        if dy>=max_pgdist|dy<=min_pgdist
            prob_map_pg(j,i)=0;
        else
            prob_map_pg(j,i)=interp1(pg_pdf(:,1), pg_pdf(:,2),dy);
        end
    end
end

if plot_flag==1
    figure;
    plot(pg_pdf(:,1),pg_pdf(:,2))
    title('pg pdf')

    figure;
    plot(mg_pdf(:,1),mg_pdf(:,2))
    title('mg pdf')

    figure;
    imagesc(g.x,g.y,prob_map_mg,[0 .03]);
    title('Migration Probability Map');
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    axis xy;
    axis equal

    figure;
    imagesc(g.x,g.y,prob_map_pg,[0 .03]);
    title('Progradation Probability Map');
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    axis xy;
    axis equal
end
      
      
