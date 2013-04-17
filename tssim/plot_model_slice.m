%plots slice number n
%code: 1: X-Z section 
%      2: Y-Z section
%by Alejandro Leiva
function plot_model_slice(surface_cube,tracklobe,nslice,code,legend_flag)
x_max = max(max(max(surface_cube)));
fgu = 1;
lobe = 1;
iter = 0;
[m n l]=size(surface_cube);
plots = find(tracklobe==-2);
if numel(plots)==0
    plots=l;
else
    plots = plots(1)-1;
end
hold;
if code==1
    title(strcat('Cross Section X-Z, Slice: ',num2str(nslice)));
    plot(surface_cube(nslice,:,1),'--b','LineWidth',3);
else
    title(strcat('Cross Section Y-Z, Slice: ',num2str(nslice)));
    plot(surface_cube(:,nslice,1),'--b','LineWidth',3);
end
legend_string(1) = {'base topo'};
%lobes
for i=2:plots(1)
    if tracklobe(i)>0
        legend_string(lobe + 1) = {strcat('lobe',num2str(tracklobe(i)))};
        lobe = lobe + 1;
        if code==1
            plot(surface_cube(nslice,:,i),'Color',[0 1 0],'LineWidth',2);        
        else
            plot(surface_cube(:,nslice,i),'Color',[0 1 0],'LineWidth',2);
        end
    iter = lobe;    
    end    
end
%fgu
for i=2:plots
    if tracklobe(i)==0
        legend_string(iter + fgu) = {strcat('fgu',num2str(fgu))};
        fgu = fgu + 1;
        if code==1
            plot(surface_cube(nslice,:,i),'Color',[1 0 0],'LineWidth',2)
        else
            plot(surface_cube(:,nslice,i),'Color',[1 0 0],'LineWidth',2)
        end  
    end
end
if legend_flag == true
    legend(legend_string);
end
if code==1
    axis([0  n  0  x_max])
    xlabel('X coord.');
    ylabel('Z coord.');
else
    axis([0  m  0  x_max])
    xlabel('Y coord.');
    ylabel('Z coord.');
end

grid on;box on;
set(gca,'LineStyle','--','PlotBoxAspectRatio',[3.5 1 1]);
set(gcf,'Color',[1 1 1]);
hold off;