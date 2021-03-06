function [np, re] = lagrange_interpolation(f,xx)

dd = divided_difference(f, xx);
np = nested_polynomial(dd, xx(1:end-1));

y = subs(f,xx);
y_lagrange = subs(np,xx);

re = double(abs(y - y_lagrange) ./ y);

fplot(f, [xx(1) xx(end)])
hold on
plot(xx, y_lagrange, 'bd')
plot(xx, y_lagrange)
hold off


end
