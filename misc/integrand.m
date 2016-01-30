function [f] = integrand(x,y,l,m,K)
    f = zeros(size(x));
    for i=1:size(x,1)
        for j=1:size(x,2)
            f(i,j) = (x(i,j)-l(1))*(y(i,j)-l(2))*mvnpdf([x(i,j),y(i,j)],m,K);
        end
    end    
end 
