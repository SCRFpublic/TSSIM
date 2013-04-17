%script to compute some statistics of surface cube output from tssim

[model_facies model_sec_var model_facies_geoeas model_sec_var_geoeas]= ...
    get_model(surface_cube,tracklobe,150,100,false,'model_output.txt',0.0,3,0.3);
model1 = 0*(model_facies<3)+1*(model_facies==3);
[L,NUM] = bwlabeln(model1,26);
geo_bodies = zeros(NUM,1);
for i=1:NUM
    geo_bodies(i) = length(find(L==i));
end
n_geo = find(geo_bodies>20);
number_geo = length(n_geo); %number of independent geobodies
for i = 1:number_geo
    disp(strcat('geobody ',num2str(i),' has ',num2str(geo_bodies(n_geo(i))),' gridblocks')); 
end

sand_blocks =sum(sum(sum(model_facies==3)))
shale_blocks = sum(sum(sum(model_facies==0)))
total_blocks =  sand_blocks + shale_blocks
sand_shale_ratio = sand_blocks/shale_blocks 
sand_percent = sand_blocks*100/total_blocks
shale_percent = shale_blocks*100/total_blocks

ninety_prc = length(find(model_sec_var_geoeas<=0.9 & model_sec_var_geoeas>0));
fifty_prc = length(find(model_sec_var_geoeas<=0.5 & model_sec_var_geoeas>0));
ten_prc = length(find(model_sec_var_geoeas<=0.1 & model_sec_var_geoeas>0));

ninety = ninety_prc*100/sand_blocks
fifty = fifty_prc*100/sand_blocks
ten = ten_prc*100/sand_blocks
