function np = lagrange_interpolation(f,xx)

dd = divided_difference(f, xx);
np = nested_polynomial(dd, xx(1:end-1));

end
