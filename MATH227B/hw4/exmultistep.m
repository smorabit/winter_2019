% Program 6.7 Multistep method
% Inputs: time interval inter,
% ic=[y0] initial condition, number of steps n,
% s=number of (multi)steps, e.g. 2 for 2-step method
% Output: time steps t, solution y
% Calls a multistep method such as ab2step.m
% Example usage: [t,y]=exmultistep([0,1],1,20,2)
function [t,y]=exmultistep(inter,ic,n,s)

    h=(inter(2)-inter(1))/n;
    % Start-up phase
    y(1,:)=ic;t(1)=inter(1);
    
    for i=1:s-1 % start-up phase, using one-step method
        t(i+1)=t(i)+h;
        y(i+1,:)=trapstep(t(i),y(i,:),h);
        f(i,:)=ydot(t(i),y(i,:));
    end
    
    for i=s:n % multistep method loop
        t(i+1)=t(i)+h;
        f(i,:)=ydot(t(i),y(i,:));
        y(i+1,:)=ab2step(t(i),i,y,f,h);
        %y(i+1,:)=am2step(t(i),i,y,f,h);
    end
%     plot(t,y, 'o');
%     xlabel("t");
%     ylabel("y(t)");
%     grid on;
end
