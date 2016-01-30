
function[y] = llvm(x)

load llvm

z = round(x(:))';

found = 0;
for i=1:1024
    if z == data(i,1:10)
        perf = data(i,11)/18;   %5.5 to 7.75
        negmem = 1-(data(i,12)-10)/20;
        found = 1;
        break
    end
end

if found == 0
    perf=0;
    negmem=0;
end
y = [perf, negmem]';


