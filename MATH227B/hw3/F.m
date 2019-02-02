% function F
function y = F(x)
x1 = x(1);
x2 = x(2);
y = zeros(2,1);
y(1) = x1^2-x2^2+2*x2; % f1(x1,x2)
y(2) = 2*x1+x2^2-6; % f2(x1,x2);
end