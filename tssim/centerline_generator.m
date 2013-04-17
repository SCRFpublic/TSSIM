% This function generates centerline connecting 'head' and 'tail'
% The centerline is a straight line
% Written by: Hongmei Li
% Date: 2007 summer

% TODO-need to used flwxiblw function to generate curvelinear curve for
% centerline which connect 'head' and 'tail'

%Input
%theta-channel orientation (assume channel is straight)
%head- the channel start cell (x,y)
%tail- channel end cell(x,y)
%g-structure containing calculated grid information

%Output
%d-structure containing
% xi - x coord interpolated along channel centerline
% yi - y coord interpolated along channel centerline
% ii - cell number in x direction corresponding with each xi
% jj - cell number in y direction corresponding with each yi
% si - incremental total distance at each cell along the centerline (last
    % entry equals total line length
   
function d=centerline_generator(head,tail,g)
N=abs(round(head(2)-tail(2)));  %HM: not sure what N is here. Head and tail(1,2) are the y-coords I think (maybe grid cell numbers actually). I think head and tail are a 1x2 matrices: [x,y].
d.xi=zeros(1,N);    % HM: in this case, d is the centerline. I think d is a set of 1xN vectors of x and y grid coords where the centerline is.
d.yi=zeros(1,N);
if head(1,1)==tail(1,1)   % straight line. HM: if the x coord of the head and tail are equal, it's just a straight line - angle =0??
   for i=1:N
       if (head(1,2)-i)>0 % just in case y gets too close to the end.
     d.xi(i)=g.x(head(1,1));    % HM: all have the x-coord of the head
     d.yi(i)=g.y(head(1,2)-i+1);   %HM: the y-coords are the y-coords of the tail minus i, which counts from 1 to the number of coords in the centerline. 
     d.ii(i)=head(1,1);     %HM: I guess g.x and g.y (thus d.xi and d.yi) are the actual spatial locations, while d.ii and d.jj are the cell numbers? 
     d.jj(i)=head(1,2)-i+1;
       end
   end
elseif head(1,2)==tail(1,2) % HM: if the y-coords are equal
    N=abs(round(head(1,1)-tail(1,1)));  %HM: N is the distance in cells from the x-loc of the head and tail.
    for i=1:N       %HM: analagous to above but flip x and y and add instead of subtract.
      if (head(1,2)-i)>0 % just in case y gets too close to the end.
     d.yi(i)=g.y(head(1,2));
     d.xi(i)=g.x(head(1,1)+i-1);
     d.jj(i)=head(1,2);
     d.ii(i)=head(1,1)+i-1;
      end
   end
else
    %draw exponential line  HM: not sure what that means.
    disty=g.y(head(1,2))-g.y(tail(1,2));    %HM: the actual distance in the y-direction between the head and tail   
    distx=g.x(tail(1,1))-g.x(head(1,1));    %HM: the actual distance in the x-direction between the head and tail
    %coeff=log((disty))/(N*g.dx);
    coeff=distx/disty;
    d.yi(1)=g.y(head(1,2));
    d.xi(1)=g.x(head(1,1));
    d.jj(1)=head(1,2);
    d.ii(1)=head(1,1);
    for i=2:N       % HM: this just creates some line from an equation. In fact, it might just be a straight line. It looks like she tried to do an exp shape but switched back to a line.
        if (head(1,2)-i)>0 % just in case y gets too close to the end.
%         if distx>0  
%            d.xi(i)=(exp((i-1)*g.dx*coeff)+g.x(head(1,1)));
%         else
%             d.xi(i)=-1*(exp((i-1)*g.dx*coeff))+g.x(head(1,1));
%         end
       d.yi(i)=g.y(head(1,2)-i+1);
       d.xi(i)=coeff*(d.yi(1)-d.yi(i))+d.xi(1);
       d.jj(i)=head(1,2)-i+1;
       d.ii(i)=round(d.xi(i)/g.dx);
        end
    end
end
d.si(1)=0;      %HM: I think this is just the line starting from the origin? Or maybe just the incremental distance along the line?
for i=2:N
    d.si(i)=d.si(i-1)+sqrt((d.xi(i)-d.xi(i-1))^2+(d.yi(i)-d.yi(i-1))^2);
end    

% figure(gcf+1);clf;
% plot(round(d.xi),round(d.yi));  
% xlim([0,g.nx*g.dx]);
% ylim([0,g.ny*g.dy]);