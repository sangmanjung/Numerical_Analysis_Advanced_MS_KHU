%% MATH7003-00: Assignment #3-(1), 2019310290 Sangman Jung _ revision 2.
clear,clc

% Initial parameters, initial conditions
h = [0.25 0.5]; % step size
Y = @(x) 1./(1+x)'; % exact solution of y'
iter = 100; % iteration number for the stop criterion

% Trapezoidal method iteration to solve y'=-y^2
% Using predictor-corrector method, we can obtain the scheme as follows.
for h_iter = 1:2 % the loop of h and 2h
    x = 0:h(h_iter):5; % apply a different step size
    y0 = 1; % initial condition
    y(1) = y0; % initialize
    for j = 1:length(x)-1 % Trapezoidal method loop
        y1 = y0-h(h_iter)*y0^2; % first term using Euler's method
        y2 = y0-2*h(h_iter)*y1^2; % second term using midpoint method
        y_p = y2; % the predictor
        for i = 1:iter % order of convergence criterion
            y_c = y0-h(h_iter)*(y0^2+y_p^2)/2; % the corrector (In case, y_c == numerical y) 
            if abs(y_c-y_p) <= 10^-7 % I will present the results per 10^-6.
                break;
            end
            y_p = y_c; % update the predictor until satisfying the criterion
        end
        y0 = y_p; % set y^(j)_{n+1} term
        y(j+1) = y_c; % update our results
    end
    fval(h_iter) = {y'}; % save the values of y have different size for y_h and y_2h
    xval(h_iter) = {x'}; % save the values of x have different size for y_h and y_2h
end

% Print out setting
y_h = cell2mat(fval(1)); % numerical y for h = 0.25
y_2h = cell2mat(fval(2)); % numerical y for h = 0.5
y_h = y_h(5:4:end);
y_2h = y_2h(3:2:end);
t = 1:5;

% Print the Table 6.8
fprintf("Table 6.8  Trapezoidal method and Richardson error estimation\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
fprintf("|   x   |   y_{2h)(x)   |   Y(x) - y_{2h}(x)   |   y_{h}(x)   |   Y(x) - y_{h}(x)   |   [y_{h}(x) - y_{2h}(x)]/3   |\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
for i = 1:length(t)
    fprintf('   %1.1f      %1.6f           %1.6f            %1.6f           %1.6f                  %1.6f\n',...
        abs([t(i) y_2h(i) Y(t(i))-y_2h(i) y_h(i) Y(t(i))-y_h(i) (y_h(i)-y_2h(i))/3]));
end
fprintf("--------------------------------------------------------------------------------------------------------------------\n");

