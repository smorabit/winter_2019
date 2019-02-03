function x = fig_newton(F, max_iter, acc)

    syms x1 x2;
    x = zeros(2,1);
    J = jacobian(F, [x1,x2]);
    i = 1;
    errors = [];
    while i < max_iter
        Ji = double(subs(J, {x1, x2}, {x(1), x(2)}));
        Fi = double(-1*subs(F, {x1, x2}, {x(1), x(2)})).';
        y = Ji \ Fi;
        x = x+y;
        err = norm(y,2);
        errors = [errors err];
        if err < acc  break; end
        i = i + 1;
    end
    if i >= max_iter
        disp("max iterations reached")
    end
    
    semilogy(errors)
end