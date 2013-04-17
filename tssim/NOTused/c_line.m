
function d=c_line(base_topo,center,sp,sae,g) %sae: sim area edge.(i,j)
% used by tssim. generate centerline connecting 'head' and 'tail'

%Input
%theta-channel orientation (assume channel is straight)
%center- the channel start cell (x,y)
%sp- channel end cell(x,y)
%sae - further point from source point on the simulation area boundary
%       sae: sim area edge.(i,j)
%g-structure containing calculated grid information

%Output
%d-structure containing
% xi - x coord interpolated along channel centerline
% yi - y coord interpolated along channel centerline
% ii - cell number in x direction corresponding with each xi
% jj - cell number in y direction corresponding with each yi
% si - incremental total distance at each cell along the centerline (last
    % entry equals total line length
    
%fix the case sp(1)<3


if sae(1)> sp(1)+ 10 
    if sae(2)==sp(2)
        sae(2) = sae(2) + 1;
        if center(2)==sp(2)
            center(2) = center(2) + 1;
        end
    end
    mid_point = ceil((center + sp)/2);
    dx = abs(center(2) - sp(2))/2;
    if dx>6
        dx = 6;
    end
    startj =mid_point(1)-dx;
    if startj < 1
        startj = 1;
    end
    min_h = min(base_topo(mid_point(1),startj:mid_point(1)+dx));
    startj =mid_point(1)-dx;
    if startj < 1
        startj = 1;
    end
    min_h = find(base_topo(mid_point(1),startj:mid_point(1)+dx)==min_h);
    mid_point(3) = mid_point(2) - dx + min_h-1;


    yi = [center(2) mid_point(2) sp(2) sae(2)];
    xi = [center(1) mid_point(1) sp(1) sae(1)];

    y2 = [center(2) mid_point(3) sp(2) sae(2)];
    x2 = [center(1) mid_point(1) sp(1) sae(1)];

    pp = interp1(y2,x2,'cubic','pp');
    xx = center(1):sae(1)/300:sae(1);
    yy = ppval(pp,xx);

    hold on
    yy = spline(x2,y2,xx);
    plot(yi,xi,'bo');
    plot(y2,x2,'go');
    plot(yy,xx);

else %not to cosider sae
    if center(2)==sp(2)
        center(2) = center(2) + 1;
    end
    mid_point = ceil((center + sp)/2);
    dx = abs(center(2) - sp(2))/2;
    if dx>6
        dx = 6;
    end
    startj =mid_point(1)-dx;
    if startj < 1
        startj = 1;
    end
    min_h = min(base_topo(mid_point(1),startj:mid_point(1)+dx));
    startj =mid_point(1)-dx;
    if startj < 1
        startj = 1;
    end
    min_h = find(base_topo(mid_point(1),startj:mid_point(1)+dx)==min_h);
    mid_point(3) = mid_point(2) - dx + min_h-1;

    yi = [center(2) mid_point(2) sp(2)];
    xi = [center(1) mid_point(1) sp(1)];

    y2 = [center(2) mid_point(3) sp(2)];
    x2 = [center(1) mid_point(1) sp(1)];

    pp = interp1(y2,x2,'cubic','pp');
    xx = center(1):sp(1)/300:sp(1);
    yy = ppval(pp,xx);

    hold on
    yy = spline(x2,y2,xx);
    plot(yi,xi,'bo');
    plot(y2,x2,'go');
    plot(yy,xx);
end    



N=abs(round(head(2)-tail(2))); 
d.xi=zeros(1,N);    
d.yi=zeros(1,N);
if head(1,1)==tail(1,1)   
   for i=1:N
       if (head(1,2)-i)>0 
     d.xi(i)=g.x(head(1,1));    
     d.yi(i)=g.y(head(1,2)-i+1);   
     d.ii(i)=head(1,1);     
     d.jj(i)=head(1,2)-i+1;
       end
   end
elseif head(1,2)==tail(1,2) 
    N=abs(round(head(1,1)-tail(1,1)));  
    for i=1:N      
      if (head(1,2)-i)>0 
     d.yi(i)=g.y(head(1,2));
     d.xi(i)=g.x(head(1,1)+i-1);
     d.jj(i)=head(1,2);
     d.ii(i)=head(1,1)+i-1;
      end
   end
else
    
    disty=g.y(head(1,2))-g.y(tail(1,2));     
    distx=g.x(tail(1,1))-g.x(head(1,1));    
    coeff=distx/disty;
    d.yi(1)=g.y(head(1,2));
    d.xi(1)=g.x(head(1,1));
    d.jj(1)=head(1,2);
    d.ii(1)=head(1,1);
    for i=2:N       
        if (head(1,2)-i)>0 
            d.yi(i)=g.y(head(1,2)-i+1);
            d.xi(i)=coeff*(d.yi(1)-d.yi(i))+d.xi(1);
            d.jj(i)=head(1,2)-i+1;
            d.ii(i)=round(d.xi(i)/g.dx);
        end
    end
end
d.si(1)=0;      
for i=2:N
    d.si(i)=d.si(i-1)+sqrt((d.xi(i)-d.xi(i-1))^2+(d.yi(i)-d.yi(i-1))^2);
end  








