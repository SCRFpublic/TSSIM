% Functions that generates lobate part of the geobody. 
% Inputs: 
%     start_point: anchor point.
%     center_point: sediment source location.
%     angles: valid angle range obtained with get_angle
%     g: grid dimensions.
%     verbose: on the fly commenting flag.
%     Input_path: path (get_dirname).
% outputs:
%     lobe_ind_map: binary 2D lobe image.
%written by H. Michaels.
%modified by Alejandro Leiva, June '09
function [lobe_ind_map]=lobe_generator(start_point,center_point,angles, g,verbose,Input_path)  
%---------------------------------------------------------------
% Get all of the CDF data from the process model
lobe_length_cdf_data=[Input_path '\CDFS\lobe_length_cdf.txt'];
lobe_width_cdf_data=[Input_path '\CDFS\lobe_width_cdf.txt'];

lobe_length_cdf=read_cdf_data(lobe_length_cdf_data,verbose);
width_cdf=read_cdf_data(lobe_width_cdf_data,verbose);

L=interp1(lobe_length_cdf(:,2), lobe_length_cdf(:,1),rand); 
W=interp1(width_cdf(:,2), width_cdf(:,1),rand);

phi=angles(ceil(rand*size(angles,1)))*pi/180;

theta=-pi/4:0.0001:pi/4;

a=L;
b=W;

if start_point(1)<center_point(1)
    xstart=(start_point(1)+1)*g.dx;
elseif start_point(1)>center_point(1)
    xstart=(start_point(1)-1)*g.dx;
else
    xstart=(start_point(1))*g.dx;
end

if start_point(2)<(g.ny-4)
    ystart=(start_point(2)+1)*g.dy;
else
    ystart=(start_point(2)+(g.ny-start_point(2)))*g.dy;
end

r=a*cos(2*theta);
x = r.*cos(theta);
max_y = max(r.*sin(theta));
y=b/2*r.*sin(theta)/max_y;
x1 = x*cos(phi)-y*sin(phi)+xstart;
y1 = x*sin(phi)+y*cos(phi)+ystart;
%lobe_thick_map=zeros(g.ny,g.nx);
lobe_ind_map=zeros(g.ny, g.nx);
x=x1;
y=y1;
% rotation_m=[cos(phi) -sin(phi);sin(phi) cos(phi)]^-1;
for j=1:g.ny
    ycurr=j*g.dy;
    x1=99999;
    x2=99999;
    ydiffvect=abs(y-ycurr);
    y_diff_L2=ydiffvect<=g.dy/2;  
    min_indices=find(y_diff_L2);
    xvect=x(min_indices);
    if (xvect)
        x1=min(xvect);
        x2=max(xvect);
    end
    for i=1:g.nx
        xcurr=i*g.dx;
        if (xcurr>=x1)&&(xcurr<=x2)
            lobe_ind_map(j,i)=1;
%             x_r=xcurr-xstart;
%             y_r=ycurr-ystart;              
%             m=rotation_m*[x_r;y_r];
%             x_r=m(1);
%             y_r=2/b*m(2)*max_y;
%             r_r=sqrt(x_r^2+y_r^2);
%             ang=atan(y_r/(x_r+0.000001));
%             proportion=r_r/(a*cos(2*ang));
%             if proportion<0
%                 proportion=0;
%             elseif proportion>1
%                 proportion=1;
%             elseif  isnan(proportion)
%                 proportion=1;
%             end
%         lobe_thick_map(j,i)=lobe_max_thick*(1-proportion);
        end
    end 
end 
   



       
     