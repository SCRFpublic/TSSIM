%distances in A-B from B
%B:sim_layer
%A:sim_area
function result=distanceAB(A,B)

binary = 0.*(A==0)+(B>0);
binary = bwdist(binary);
binary = 1-(1/max(max(binary))).*binary;
result = binary;
result = result.*(A>0);
x = find(A>0);
j=min(result(x));
result = (result-j).*(result>0/(1-j));
result = result/max(max(result));

