function x = NewtonMethod(funcF,JacobianF,n)
    
F = funcF;
J = JacobianF;
x = zeros(n,1); % set initial to (0,...,0)
Iter = 1;
MaxIter = 100;
TOL = 1e-5;

while Iter < MaxIter
    y = J(x)\(-F(x));
    x = x+y;
    if norm(y,2)<TOL break; end
end
if Iter >= MaxIter
    disp('Maximum number of iteration exceeded!');
end
end
