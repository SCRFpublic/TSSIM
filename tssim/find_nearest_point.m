% Find nearest points on centerline defined in d to every point on grid 
% defined by g
% d - structure containing
% xi - x coord interpolated along channel centerline
% yi - y coord interpolated along channel centerline
% si - along channel centerline length
% g - structure containing calculated grid information
% verbose - 1:comments otherwise no comments
%
% Return
% across_channel_length - least squares distance from each point on
% centerline to every point on grid
% index - index at every point on grid of nearest point on centerline
% Written by Jonathan Steward

function [across_channel_distance,index]=find_nearest_point(d,g,verbose)



% Length of data
N=length(d.xi);

% if(verbose)
%     fprintf('find_nearest_point: Finding nearest point on centerline to every point on grid ');
% end
% tic;

% Set up arrays
tmp_dist=zeros(1,N);
across_channel_distance=zeros(g.ny,g.nx);
index=zeros(g.ny,g.nx);

for i=1:g.nx 
    for j=1:g.ny
        % add a random component to make sure distances are different
        tmp_dist=0.01*rand(1)+(((d.xi-g.x(i)).^2)+((d.yi-g.y(j)).^2)).^0.5;
        mtd=min(tmp_dist);
        %temp_dist=(((d.xi-g.x(i)).^2)+((d.yi-g.y(j)).^2)).^0.5;
        %temp_dist =((((d.ii-i).^2+(d.jj-j).^2).^0.5))*g.dx;
        ii=find(tmp_dist==mtd);
        across_channel_distance(j,i)=tmp_dist(ii(1));
        index(j,i)=ii(1);
    end
end

% if(verbose==1)
%     fprintf('in %.2f seconds\n',toc);
%     figure(gcf+1);clf;
%     imagesc(g.x,g.y,across_channel_distance);
%     colorbar;
%     hold on;
%     plot(round(d.xi),round(d.yi),'w');
%     hold off;
%     xlabel('x [m]');
%     ylabel('y [m]');
%     axis xy;
%     title('find nearest point');
% end
