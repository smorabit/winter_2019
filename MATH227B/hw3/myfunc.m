function F = myfunc(x)

aMin = 3;
aMax = 10;
aDeg = 3;
bMin = 2;
bMax = 12;
bDeg = 5;
eX = 9;
eY = 11;
n = 1;

F(1) = aMin - x(1)*aDeg + (aMax - aMin) * (x(2)^n)/(eY^n + x(2)^n);
F(2) = bMin - x(2)*bDeg + (bMax - bMin) * (x(1)^n)/(eX^n + x(1)^n);

end
