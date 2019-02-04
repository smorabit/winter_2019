function x = math(start)

syms x1 x2;

aMin = 0.1;
aMax = 8;
aDeg = 5;
bMin = 0.1;
bMax = 10;
bDeg = 5;
eX = 0.5;
eY = 0.5;
n = 2;

f1 = aMin - x1*aDeg + (aMax - aMin) * (x2^n)/(eY^n + x2^n);
f2 = bMin - x2*bDeg + (bMax - bMin) * (x1^n)/(eX^n + x1^n);
F = [f1, f2];

x = fig_newton(F,100,1e-5,start); 
J = jacobian(F, [x1,x2]);

%evaluate the jacobian at the fixed point to determine stability:
J_fp = double(subs(J, {x1,x2}, x.'));
Trace = trace(J_fp)
Det = det(J_fp)


end
