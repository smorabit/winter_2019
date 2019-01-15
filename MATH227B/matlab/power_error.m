function o = power_error(f, x, desired_error)
    
    %note: x parameter should be between [0,1]

    o = 0;
    error = realmax;
    
    while error > desired_error
        o = o + 1;
        
        %compute taylor polynomial
        T = taylor(f, "Order", o);
        
        %compute error between taylor and real polynomial
        error = relative_error(double(subs(f,x)), double(subs(T,x)));
end