% This function generate channel thickness map using a parabola function
% Written by: Hongmei Li
% Date: 2007 Summer
% Modified by Holly Michael to use both the erosional and depositional
% thicknesses, putting them on the same thickness map, and making it change
% along the length of the channel.

%across_chan_dist- the distance map relative to channel center line
%chan_detrend_totalL_to_width_cdf - cdf of total channel length to width,
    %trend is 0.0249x + 0.6064
%chan_h- channel max. thickness

function [chan_thick_map,chan_ind_map]=channel_shape_generator(across_chan_dist,chan_detrend_totalL_to_width_cdf,chan_h,d,g)
% 
chan_length=d.si(end);
chan_detrend_LtoW=interp1(chan_detrend_totalL_to_width_cdf(:,2), chan_detrend_totalL_to_width_cdf(:,1),rand);
chan_LtoW=(0.0249*chan_length+.6064)+chan_detrend_LtoW;
chan_w=chan_length/chan_LtoW;

% Let erosional length and width be a random proportion of the length and
    % width: U (0.5,1)
  %  upper=1;
  %  lower=0.5;
% e_chan_L=chan_length;   % make them equal so there's no space between channel and lobe
% e_chan_W=(lower+(upper-lower)*rand)*chan_w;

chan_thick_map=zeros(g.ny,g.nx);
%chan_eros_map=zeros(g.ny,g.nx);
chan_ind_map=zeros(g.ny,g.nx);
% parabola
for i=d.jj(1):-1:d.jj(end)
    for j=1:g.nx
        %along_chan_dist=sqrt((d.xi(1)-g.x(j))^2+(d.yi(1)-g.y(i))^2);
        %along_chan_dist=sqrt((d.xi(g.ny-i)-g.x(j))^2+(d.yi(g.ny-i)-g.y(i))^2);
      if across_chan_dist(i,j)<=chan_w/2 %&& along_chan_dist<chan_length 
     
          %chan_thick_map(i,j)=-1*chan_h*(1-(2*across_chan_dist(i,j)/chan_w)^2);
         % chan_thick_map(i,j)=-1*chan_e*(1-(2*across_chan_dist(i,j)/chan_w)^2)*(1-along_chan_dist/d.si(end))+chan_h*(1-(2*across_chan_dist(i,j)/chan_w)^2);
         
         % Make the channel have no erosion for now.
        
         %across_prop(find(across_prop>1.0))=1.0;
         %chan_thick_map(i,j)=((acros_prop)^2)*chan_h*(1-(2*across_chan_dist(i,j)/chan_w)^2);
         chan_thick_map(i,j) = chan_h*sqrt(1-across_chan_dist(i,j)^2/(chan_w/2)^2);
          % indicate where the channels are, and indicate the location of
          % erosion. Chose -2 so that if added to the lobe map, will still
          % be negative.
          if chan_thick_map(i,j)>0
              chan_ind_map(i,j)=1;
          elseif chan_thick_map(i,j)<0
              chan_ind_map(i,j)=-2;
          end
      end
%       if across_chan_dist(i,j)<=e_chan_W/2 && along_chan_dist<e_chan_L 
%           % just make erosional thickness constant in the channel
%           chan_eros_map(i,j)=chan_e;
%       end
   end
end

% figure;
% imagesc(g.x,g.y,chan_eros_map);
% title('channel erosion map');
% colorbar;
% axis xy;
