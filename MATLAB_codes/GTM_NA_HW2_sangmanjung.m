%% MATH7003-00: Assignment #2, 2019310290 Sangman Jung
clear,clc

% parameters, initial conditions
h = 0.25; % step size
x = 0:0.25:2.25;
lambda = -1; % parameter. the case of lambda < 0
y(1) = 1; % initial condition
y(2) = y(1)+h*lambda*y(1); % obtained using Euler's method

% general solution of midpoint scheme
r_0 = h*lambda + sqrt(1+h^2*lambda^2); % root r0
r_1 = h*lambda - sqrt(1+h^2*lambda^2); % root r1
beta_0 = (y(2)-r_1*y(1))/(r_0-r_1); % coefficient beta0
beta_1 = (y(1)*r_0-y(2))/(r_0-r_1);% coefficient beta1

% Midpoint method for the equation { y' = lambda*y , y(0) = 1 }
fprintf("Table 6.6  Example 1 of midpoint method instability\n");
fprintf("----------------------------------------------------\n");
fprintf("|   x_{n}   |   y_{n}   |   Y(x_{n})   |   Error   |\n");
fprintf("----------------------------------------------------\n");
for n = 2:length(x)
    y(n+1) = y(n-1)+2*h*lambda*y(n); % midpoint method
    Y(n) = exp(-x(n)); % exact solution of y'
    Error(n) = Y(n)-y(n);
    fprintf('     %1.2f       %1.4f       %1.4f       %1.4f    \n',...
        [x(n) y(n) Y(n) Error(n)]);
end
fprintf("----------------------------------------------------\n");

% Parasitic solution of the midpoint method
fprintf('\n  *** Parasitic solution of the midpoint method ***\n');
fprintf("----------------------------------------------------\n");
for n = 1:length(x)
    para_sol(n) = beta_1*r_1^n;
    fprintf('  x_{n} = %1.2f  and  beta_{1}*r_{1}^(%1d) = %1.4f\n',[x(n) n para_sol(n)]);
end
fprintf("----------------------------------------------------\n");

% Graphs
subplot(1,2,1)
plot(x,y(1:10),'g-o',x,Y,'k-o','LineWidth',1.5);
xlim([0 x(end)]);ylim([-0.2 1.2]);
xlabel('x');
legend('Numerical','Exact');
title('Numerical results vs Exact results');
tt = get(gca,'title');
tt.FontWeight = 'bold';
% grid on
subplot(1,2,2)
plot(x,para_sol,'r-o',x,Error,'b-o','LineWidth',1.5);
xlim([0 x(end)]);ylim([-0.2 0.2]);
xlabel('x');
legend('Parasitic','Error');
title('Parasitic solution vs Error');
tt = get(gca,'title');
tt.FontWeight = 'bold';
% grid on