%% MATH7003-00: Assignment #1, 2019310290 Sangman Jung
clear,clc

% parameters, initial conditions
K = ellipke(1/2); % the complete elliptical integral of the first kind
h = [0.2 0.1]; % step size
y1(1) = pi/2; % initial value of y1
y2(1) = 0; % initial value of y2

% Euler method for the pendulum equation
fprintf("Table 6.5  Euler's method for example (6.2.57)\n");
fprintf("-------------------------------------------------------------------------------\n");
fprintf("| h | x_{n} | y_{1,n} | Y_{1}(x_{n}) | Error | y_{2,n} | Y_{2}(x_{n}) | Error |\n");
fprintf("-------------------------------------------------------------------------------\n");
for h_iter = 1:2
    x = 0:h(h_iter):1;
    for n = 1:length(x)
        y1(n+1) = y1(n)+h(h_iter)*y2(n);
        y2(n+1) = y2(n)-h(h_iter)*sin(y1(n));
        [SN(n),CN(n),DN(n)] = ellipj(K-x(n),1/2); % the Jacobi elliptic functions
        Y1(n) = 2*asin(sqrt(2)/2*SN(n)); % the exact solution of y1
        Y2(n) = -2*CN(n)*DN(n)/sqrt(2-SN(n)^2); % the exact solution of y2
        Error1(n) = Y1(n)-y1(n);
        Error2(n) = Y2(n)-y2(n);
        fprintf('%1.1f   %1.1f    %1.4f       %1.4f    %1.4f   %1.5f   %1.6f   %1.6f\n',...
            [h(h_iter) x(n) y1(n) Y1(n) Error1(n) y2(n) Y2(n) Error2(n)]);
    end
    fprintf("-------------------------------------------------------------------------------\n");
end