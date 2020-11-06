%% MATH7003-00: Assignment #7, 2019310290 Sangman Jung.
clear, close, clc
% functions
f = @(x,y) (1./(1+x.^2)) - 2*y.^2; % ODE
Y = @(x) x./(1+x.^2); % exact solution
ktable = 3:7;
K = 1./2.^(ktable); % additional topics in HW#7
h = [0.25 0.5 K 2*K]; % step size h (1x12)

% Runge-Kutta method iteration
for h_iter = 1:length(h) % h = 0.25, 0.5 & h = (1/2)^k
    x = 0:h(h_iter):10; % variable x per step size
    y = zeros(1,length(x)); % allocate the numerical solution y
    y(1) = 0; % initial condition
    for n = 1:length(x)-1 % Runge-Kutta iteration
        V1(n) = f(x(n),y(n));
        V2(n) = f(x(n)+h(h_iter)/2,y(n)+(h(h_iter)*V1(n))/2);
        V3(n) = f(x(n)+h(h_iter)/2,y(n)+h(h_iter)/2*V2(n));
        V4(n) = f(x(n)+h(h_iter),y(n)+h(h_iter)*V3(n));
        y(n+1) = y(n) + h(h_iter)*(V1(n)+2*V2(n)+2*V3(n)+V4(n))/6;
    end
    xvals(h_iter) = {x'}; % save the values of x for h and 2h
    yvals(h_iter) = {y'}; % save the values of y for h and 2h
end

% print out setting for Table 6.24
xh = cell2mat(xvals(1)); yh = cell2mat(yvals(1));
x2h = cell2mat(xvals(2)); y2h = cell2mat(yvals(2));

% print out setting for h = (1/2)^k
xk3 = cell2mat(xvals(3)); yk3 = cell2mat(yvals(3));
xk4 = cell2mat(xvals(4)); yk4 = cell2mat(yvals(4));
xk5 = cell2mat(xvals(5)); yk5 = cell2mat(yvals(5));
xk6 = cell2mat(xvals(6)); yk6 = cell2mat(yvals(6));
xk7 = cell2mat(xvals(7)); yk7 = cell2mat(yvals(7));
x2k3 = cell2mat(xvals(8)); y2k3 = cell2mat(yvals(8));
x2k4 = cell2mat(xvals(9)); y2k4 = cell2mat(yvals(9));
x2k5 = cell2mat(xvals(10)); y2k5 = cell2mat(yvals(10));
x2k6 = cell2mat(xvals(11)); y2k6 = cell2mat(yvals(11));
x2k7 = cell2mat(xvals(12)); y2k7 = cell2mat(yvals(12));

% find the index of Table 6.24
Table = 2:2:10;
ind_xh = zeros(1,length(Table));
ind_x2h = zeros(1,length(Table));
for i = Table
    ind_xh(i) = find(xh == i);
    ind_x2h(i) = find(x2h == i);
end
ind_xh = nonzeros(ind_xh);
ind_x2h = nonzeros(ind_x2h);

% find the index of h = (1/2)^k
kx = 2; % we want to see only the ratio at x = 2.0
ind_xk3 = find(xk3 == kx); ind_xk4 = find(xk4 == kx);
ind_xk5 = find(xk5 == kx); ind_xk6 = find(xk6 == kx);
ind_xk7 = find(xk7 == kx);
ind_x2k3 = find(x2k3 == kx); ind_x2k4 = find(x2k4 == kx);
ind_x2k5 = find(x2k5 == kx); ind_x2k6 = find(x2k6 == kx);
ind_x2k7 = find(x2k7 == kx);

% calculate the columns for Table 6.24
yh_table = yh(ind_xh);
Error_h = Y(xh(ind_xh))-yh(ind_xh);
Error_2h = Y(x2h(ind_x2h))-y2h(ind_x2h);
Ratio = Error_2h./Error_h;
Richardson = (yh(ind_xh) - y2h(ind_x2h))/15;

% calculate the columns for h = (1/2)^k
yk_t = [yk3(ind_xk3) yk4(ind_xk4) yk5(ind_xk5) yk6(ind_xk6) yk7(ind_xk7)];
y2k_t = [y2k3(ind_x2k3) y2k4(ind_x2k4) y2k5(ind_x2k5) y2k6(ind_x2k6) y2k7(ind_x2k7)];

Error_k = [Y(xk3(ind_xk3))-yk_t(1) Y(xk4(ind_xk4))-yk_t(2) ...
    Y(xk5(ind_xk5))-yk_t(3) Y(xk6(ind_xk6))-yk_t(4) ...
    Y(xk7(ind_xk7))-yk_t(5)];

Error_2k = [Y(x2k3(ind_x2k3))-y2k_t(1) Y(x2k4(ind_x2k4))-y2k_t(2) ...
    Y(x2k5(ind_x2k5))-y2k_t(3) Y(x2k6(ind_x2k6))-y2k_t(4)...
    Y(x2k7(ind_x2k7))-y2k_t(5)];

Ratio_k = [Error_2k(1)/Error_k(1) Error_2k(2)/Error_k(2) Error_2k(3)/Error_k(3) ...
    Error_2k(4)/Error_k(4) Error_2k(5)/Error_k(5)];

% Table 6.24
fprintf('Table 6.24  Example of Runge-Kutta method (6.10.21)\n');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('   x         y_{h}(x)    Y(x)-y_{h}(x)    Y(x)-y_{2h}(x)    Ratio    [y_{h}(x)-y_{2h}(x)]/15\n');
fprintf('---------------------------------------------------------------------------------------------\n');
for t = 1:length(Table)
    fprintf('  %2s       % 1.8f     % 1.1e          % 1.1e        % 2.0f            % 1.1e\n',...
        num2str(Table(t)),yh_table(t),Error_h(t),Error_2h(t),Ratio(t),Richardson(t));
end
fprintf('---------------------------------------------------------------------------------------------\n\n\n');

% The study of order of convergence at x = 2.0
fprintf('Table 6.24(HW#7) Order of convergence at x = 2.0\n');
fprintf('---------------------------------------------------------------------------------\n');
fprintf('    k    h = (1/2)^k    y_{h}(x)    Y(x)-y_{h}(x)    Y(x)-y_{2h}(x)     Ratio    \n');
fprintf('---------------------------------------------------------------------------------\n');
for T = ktable
    fprintf('    %s      %1.4f     % 1.8f     % 1.1e         % 1.1e         % 2.0f  \n',...
        num2str(T),h(T),yk_t(T-2),Error_k(T-2),Error_2k(T-2),Ratio_k(T-2));
end
fprintf('---------------------------------------------------------------------------------\n');

% Ratio graph
for i = ktable
    figtxt(i-2) = {sprintf('k=%d',i)};
end
figure('Name','Convergence analysis for R-K method')
plot(h(3:7),Ratio_k,'k-*','MarkerEdgeColor','r','MarkerSize',8); ylim([15.5 20]);
xticklabels('');
ax=gca; ax.YGrid='on';
xlabel('h = (1/2)^{k}'); ylabel('Ratio');
title('Order of convergence for R-K method at x = 2.0');
text(h(3:7),Ratio_k-0.19,figtxt,'FontSize',11)
hold on
line([0 0.14],[16 16],'Color','b','LineStyle','--','LineWidth',0.9);
legend('Numerical Ratio','Theoretical Ratio','Location','northwest');
hold off
