function[ y ] = noc( x )

width           = ceil(x(1)*7);
complexity      = round(x(2)); 
fifo            = ceil(x(1)*4);

multiplier_ind  = ceil(x(1)*7);
multiplier_list = [1,2,5,10,20,50,100];
multiplier      = multiplier_list(multiplier_ind);

load noc.mat

found = 0;
for i=1:259
    if width==inputs(i,1) && complexity==inputs(i,2) && fifo==inputs(i,3) && multiplier==inputs(i,4)
        found = 1;
        neg_eng  = 2.5-energy(i)/4;
        inv_time = (invtime(i)-4.3)*1.2;
        break
    end
end

if found == 0
    neg_eng  = 0;
    inv_time = 0;
end

y = [neg_eng, inv_time]';
