function e = lagrange_error(f, xx)

y = subs(f, xx);

%compute lagrange interpolation:
np = subs(lagrange_iterpolation(f,xx));



end