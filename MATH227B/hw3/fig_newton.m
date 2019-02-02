%F: a vector of functions
function x = fig_newton(F, max_iter, acc)

    syms x1 x2;
    x = zeros(2,1);
    J = jacobian(F, [x1,x2]);
    i = 1;
    
    while i < max_iter
        Ji = double(subs(J, {x1, x2}, {x(1), x(2)}));
        Fi = double(-1*subs(F, {x1, x2}, {x(1), x(2)})).';
        y = Ji \ Fi;
        x = x+y;
        if norm(y,2) < acc break; end
        i = i + 1;
    end
    if i >= max_iter
        disp("max iterations reached")
    end

end