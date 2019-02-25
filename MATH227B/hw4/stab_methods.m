function stab_methods( npts, method )
%
% This function makes plots of stability regions 
% of various methods using the "boundary locus method".
% For details, see:
% 
% "Finite difference methods for ODEs and PDEs" by
% R. J. LeVeque (p.g. 162).
%
% In brief, we plot the roots of the stability polynomial. 
% We assume that the complex number \zeta is represented
% by \zeta_{j} = |r_{j}| * exp(i*theta). Since |\zeta_{j}|<=1
% for stability, then \zeta_{j} = exp(i*theta). Furthermore,
%
%      z = \rho(e^{i*\theta}) / \sigma(e^{i*\theta}).
%
% To establish stability, pick a random point on \mathbb{C} 
% and check whether the root condition is satisfied. This
% area is shadowed only!
% 
%
% Inputs: 1)   npts: number of points along theta.
%         2) method: pick the method of interest.
%           
% Output: The figure corresponding to the stability region. 
%         The latter is highlighted with a shadowed gray area.
%
% Method available:
%    1) method = 'AB1'    (Forward Euler's).
%    2) method = 'BEuler' (Backward Euler's).
%    3) method = 'AB2'    (Two-step Adams-Bashforth).
%    4) method = 'AB3'    (Three-step Adams-Bashforth).
%    5) method = 'AB4'    (Four-step Adams-Bashforth).
%    6) method = 'AB5'    (Five-step Adams-Bashforth).
%    7) method = 'AM2'    (Two-step Adams-Moulton).
%    8) method = 'AM3'    (Three-step Adams-Moulton).
%    9) method = 'AM4'    (Four-step Adams-Moulton).
%
%

% Main-setup:
theta = linspace(0,2*pi,npts); % Number of points along $\theta$.
clr = [221 221 221]/255;       % User-defined color.
zeta = exp(1i*theta); 
 
figure;
set(gca,'fontsize',20,'fontname','times'); 
switch(method)
    
    case('AB1')    % Forward Euler.
        % Solution to the Characteristic Polynomial:
        z = zeta - 1; 
        plot(real(z),imag(z),'-k','linewidth',2); hold on;
        % Stable inside! Fill-in it!
        fill(real(z),imag(z),clr); 
        axis([-2 2 -2 2]);
        ht = title('Forward Euler method (AB1)');

    case('BEuler') % Backward Euler.
        % Solution to the Characteristic Polynomial:
        z = 1 - conj(zeta);
        % Fill the region with clr color!
        set(gca,'color',clr);
        hold on;
        plot(real(z),imag(z),'-k','linewidth',2); hold on;
        fill(real(z),imag(z),'w'); % Create a "white hole"! 
        axis([-2 2 -2 2]);
        ht = title('Backward Euler method');
   
    case('AB2')    % r = 2 Adams-Bashforth. 
        pi_ab2 = @(z) 2 *( z.^2 - z )./ ( 3 * z - 1 );
          vab2 = pi_ab2(zeta);
        plot(real(vab2),imag(vab2),'-k','linewidth',2); hold on;
        fill(real(vab2),imag(vab2),clr); 
        axis([-3 1 -2 2]);
        ht = title('Two-step Adams-Bashforth');
        
    case('AB3')    % r = 3 Adams-Bashforth.
        pi_ab3 = @(z) 12 *( z.^3 - z.^2 )./ ( 5 - 16 * z + 23 * z.^2 );
          vab3 = pi_ab3(zeta);
        plot(real(vab3),imag(vab3),'-k','linewidth',2); hold on;
        fill(real(vab3),imag(vab3),clr); 
        axis([-3 1 -2 2]);
        ht = title('Three-step Adams-Bashforth');
        
    case('AB4')    % r = 4 Adams-Bashforth.
        pi_ab4 = @(z) 24 *( z.^4 - z.^3 )./ ...
          ( -9 + 37 * z - 59 * z.^2 + 55 * z.^3 );
          vab4 = pi_ab4(zeta);
          plot(real(vab4),imag(vab4),'-k','linewidth',2); hold on;
          fill(real(vab4(real(vab4)<=0)),imag(vab4(real(vab4)<=0)),clr);
          axis([-3 1 -2 2]);
          ht = title('Four-step Adams-Bashforth');

    case('AB5')  % r = 5 Adams-Bashforth.
        pi_ab5 = @(z) 720 *( z.^5 - z.^4 )./ ...
          ( 251 - 1274 * z + 2616 * z.^2 - 2774 * z.^3 + 1901 * z.^4 );
          vab5 = pi_ab5(zeta);
          plot(real(vab5),imag(vab5),'-k','linewidth',2); hold on;
          fill(real(vab5(real(vab5)<=0)),imag(vab5(real(vab5)<=0)),clr);
          axis([-3 1 -2 2]);
          ht = title('Five-step Adams-Bashforth');
    
    case('AM2')    % r = 2 Adams-Moulton.  
        pi_am2 = @(z) 12 *( z.^2 - z )./ ( -1 + 8 * z + 5 * z.^2 );
          vam2 = pi_am2(zeta);
          plot(real(vam2),imag(vam2),'-k','linewidth',2); hold on;
          fill(real(vam2),imag(vam2),clr); 
          axis([-7 1 -4 4]);
          ht = title('Two-step Adams-Moulton');

    case('AM3')    % r = 3 Adams-Moulton.  
        pi_am3 = @(z) 24 *( z.^3 - z.^2 )./ ( 1 - 5 * z + 19 * z.^2 + 9 * z.^3 );
          vam3 = pi_am3(zeta);
          plot(real(vam3),imag(vam3),'-k','linewidth',2); hold on;
          fill(real(vam3),imag(vam3),clr); 
          axis([-7 1 -4 4]);
          ht = title('Three-step Adams-Moulton');
    case('AM4')    % r = 4 Adams-Moulton.
        pi_am4 = @(z) 720 *( z.^4 - z.^3 )./ ...
         ( -19 + 106 * z - 264 * z.^2 + 646 * z.^3 + 251 * z.^4);
        vam4 = pi_am4(zeta);
        plot(real(vam4),imag(vam4),'-k','linewidth',2);
        hold on;
        fill(real(vam4),imag(vam4),clr); 
        axis([-7 1 -4 4]);
        ht = title('Four-step Adams-Moulton');
end
grid on;
xlabel('$Re(z)$','interpreter','latex');
ylabel('$Im(z)$','interpreter','latex');
set(gca,'fontsize',20,'fontname','times');
set(ht,'fontsize',16,'fontname','times'); 
end