% gets contour of matrix a, 2D geobody.
% Input:
%     a: 2D matrix.
% Output: 
%     result: 2D matrix with contour of a.
%written by Alejandro D. Leiva, June '09.
function result = cont(a)
    a(a~=0)=1;
    [m n]=size(a);
    result=zeros(m,n);
    result1=zeros(m,n);
    for j=1:n
            if a(1,j)~=0
                result1(1,j)=1;
            end
            if a(m,j)~=0
                result1(m,j)=1;
            end
    end
    for i=1:m
            if a(i,1)~=0 
                result1(i,1)=1;
            end
            if a(i,n)~=0 
                result1(i,n)=1;
            end
     
    end
    for i=2:m-1
        for j=2:n-1
            if a(i,j)~=a(i,j+1) && a(i,j+1)~=0
                result(i,j)=1;
            elseif a(i,j)~=a(i,j-1) && a(i,j-1)~=0
                result(i,j)=1;
            elseif a(i,j)~=a(i+1,j) && a(i+1,j)~=0
                result(i,j)=1;
            elseif a(i,j)~=a(i-1,j) && a(i-1,j)~=0
                result(i,j)=1;
            else 
                result(i,j)=0;       
            end
        end
    end
    result=result1+result;
end

    
    
    
    

