%F: a vector of functions
function sol = fig_newton(F, max_iter, tol)
    syms x1 x2
    x = zeros(2,1);
    
    %compute jacobian of these functions
    J = jacobian(F, [x1,x2]);
    
    i = 1;
    while i < max_iter
        y = (subs(J, {x1, x2}, {x(1), x(2)}))\(-subs(F, {x1, x2}, {x(1), x(2)}));
        x = x+y;
        if norm(y,2)<tol break; end
        i = i + 1;
    end
    

end