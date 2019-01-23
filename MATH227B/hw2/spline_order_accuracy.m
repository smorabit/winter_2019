% input:
% 
% f: a function
% a: lower bound to sample function
% b: upper bound to sample function
% n: number of data points to sample
% boundary: true = use slope information for boundary

% note: y in the below function can have 


function p = spline_order_accuracy(f,a,b,n, boundary, plot)
    
    %boundary conditions: start & end slope
    start_slope = double(subs(diff(f),a));
    end_slope = double(subs(diff(f),b));
    
    %n random uniform numbers on the interval [a,b]
    x_sim = a + (b-a) * rand(n,1);
    if boundary
        y_sim = [start_slope double(subs(f,x_sim)).' end_slope].';
    else
        y_sim = double(subs(f,x_sim));
    end
    %compute spline
    h = (b - a) / (n-1); % h is step size 
    xx = a:h:b;
    y = double(subs(f,xx)); 
    yy = spline(x_sim, y_sim, xx);
    
    %compute spline for 2n:
    x_sim2 = a + (b-a) * rand(2*n,1);
    if boundary
        y_sim2 = [start_slope double(subs(f,x_sim2)).' end_slope].';
    else
        y_sim2 = double(subs(f,x_sim2));
    end
    yy2 = spline(x_sim2, y_sim2, xx);
    
    %plot
    if plot
        if boundary
            plot(x_sim,y_sim(2:end-1),'o',xx,yy);
        else
            plot(x_sim,y_sim,'o',xx,yy);
        end
        grid on;
    end
    
    
    %compute relative error between simulated data and spline:
    e = norm(abs(y - yy), Inf);
    e2 = norm(abs(y - yy2), Inf);
    
    p = log(e/e2)/log(2);
    

end

