function x = math

syms x1 x2;

aMin = 3;
aMax = 10;
aDeg = 3;
bMin = 2;
bMax = 12;
bDeg = 5;
eX = 9;
eY = 11;
n = 1;

f1 = aMin - x1*aDeg + (aMax - aMin) * (x2^n)/(eY^n + x2^n);
f2 = bMin - x2*bDeg + (bMax - bMin) * (x1^n)/(eX^n + x1^n);
F = [f1, f2];

x = fig_newton(F,100,1e-5); 
J = jacobian(F, [x1,x2]);

%evaluate the jacobian at the fixed point to determine stability:
J_fp = double(subs(J, {x1,x2}, x.'))
trace(J_fp)
det(J_fp)


end
