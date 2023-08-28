
function [Ir_mon] = sum2month(Ir)
month1 = [31,28,31,30,31,30,31,31,30,31,30,31];
month2 = [31,29,31,30,31,30,31,31,30,31,30,31];
len = size(Ir,3);
if len < 366
    month = month1;
else
    month = month2;
end
n = 1;
Ir_mon = nan(size(Ir,1),size(Ir,2),12);
for j = 1 : 12
    x = Ir(:,:,n:n+month(j)-1);
    if all(isnan(x))
        Ir_mon(:,:,j) = NaN;
    else
        Ir_mon(:,:,j) = sum(x,3,'omitnan');
    end
    n = n+month(j);
end
end