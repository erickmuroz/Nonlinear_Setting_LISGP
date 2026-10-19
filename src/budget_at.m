function n = budget_at(x, med, level)
k = find(med <= level, 1);
if isempty(k)
    n = NaN;
elseif k == 1
    n = x(1);
else
    t = (log(level) - log(med(k-1))) / (log(med(k)) - log(med(k-1)));
    n = exp(log(x(k-1)) + t*(log(x(k)) - log(x(k-1))));
end
end
