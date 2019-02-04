function F = myfunc(x)

aMin = 0.1;
aMax = 8;
aDeg = 5;
bMin = 0.1;
bMax = 10;
bDeg = 5;
eX = 0.5;
eY = 0.5;
n = 2;

F(1) = aMin - x(1)*aDeg + (aMax - aMin) * (x(2)^n)/(eY^n + x(2)^n);
F(2) = bMin - x(2)*bDeg + (bMax - bMin) * (x(1)^n)/(eX^n + x(1)^n);

end
