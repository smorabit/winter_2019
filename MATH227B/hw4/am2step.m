function z=am2step(t,i,y,f,h)
    %one step of the Adams-Moulton 2-step method
    z=y(i,:)+h*(f(i-1,:)/2-f(i,:)/2);
end
