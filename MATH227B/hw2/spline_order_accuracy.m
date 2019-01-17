% input:
% 
% f: a function
% a: lower bound to sample function
% b: upper bound to sample function
% n: number of data points to sample


function p = spline_order_accuracy(f,a,b,n)

    %n random uniform numbers on the interval [a,b]
    x_sim = a + (b-a) * rand(n,1);
    y_sim = double(subs(f,x_sim));
    
    %compute spline
    h = (b - a) / (n-1); % h is step size I think??? idk
    xx = a:h:b;
    y = double(subs(f,xx)); % actual values at points that spline is interpolating at
    yy = spline(x_sim, y_sim, xx);
    
    %compute spline for 2n:
    x_sim2 = a + (b-a) * rand(2*n,1);
    y_sim2 = double(subs(f,x_sim2));
    yy2 = spline(x_sim2, y_sim2, xx);
    
    %plot
    plot(x_sim,y_sim,'o',xx,yy)
    
    %compute relative error between simulated data and spline:
    re = abs((y - yy)./y);
    re2 = abs((y -yy2)./y);
    
    %compute p:
    p = mean(log(re ./ re2) / log(2));
    
    
    
    
    
    
    
end