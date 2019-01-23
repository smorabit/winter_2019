function dd = divided_difference(f, x)

dd = zeros(length(x), length(x));
dd(:, 1) = double(subs(f,x));

for i = 2:length(x)
    for j = 2:i
        dd(i,j) = (dd(i,j-1) - dd(i-1,j-1)) / (x(i) - x(i-j+1));
    end
end
dd = diag(dd);
end