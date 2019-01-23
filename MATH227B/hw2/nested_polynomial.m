function p_x = nested_polynomial(a, xx)

syms z

p_x = a(1);
for i = 2:length(xx)+1
    term = 1;
    for j = 1:i-1
        term = term*(z-xx(j));
    end
    term = term * a(i);
    p_x = p_x + term;
end
end