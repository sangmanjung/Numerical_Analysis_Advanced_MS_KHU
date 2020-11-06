%% MATH7003-00: Assignment #3-(1), 2019310290 Sangman Jung _ revision 0.
clear,clc

% parameters, initial conditions
h = [0.25 0.5]; % step size
y(1) = 1; % initial value of y 

% Trapezoidal method iteration to solve y'=-y^2
% Using predictor-corrector method, we can obtain the scheme as follows.
for h_iter = 1:2 % the loop of h and 2h
    x = 0:h(h_iter):5;
    y = ones(size(x)); % reset the computed results before
    y(2) = y(1)-h(h_iter)*y(1)^2; % define the term for midpoint, using Euler's method
    y_p(1)=y(1)-2*h(h_iter)*y(2)^2; % initial guess of y_p (midpoint method)
    for j =1:length(x)-1 % Trapezoidal method loop
        y_p(j+1)=y_p(j)-(y(j)*(h(h_iter)*y(j)/2-1)+y_p(j)...
            *(1+h(h_iter)*y_p(j)/2))/(1+h(h_iter)*y_p(j)); % predictor (Newton's method)
        y_c(j+1)=y(j)-h(h_iter)*(y(j)^2+y_p(j+1)^2)/2; % corrector
        y(j+1) = y_c(j+1); % update
    end
    fval(h_iter) = {y'}; % save the values have different size for y_h and y_2h
end

% print out setting
y_h = cell2mat(fval(1)); % numerical y for h = 0.25
y_2h = cell2mat(fval(2)); % numerical y for h = 0.5
y_h = y_h(5:4:end);
y_2h = y_2h(3:2:end);
t = 1:5;

% print the Table 6.8
fprintf("Table 6.8  Trapezoidal method and Richardson error estimation\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
fprintf("|   x   |   y_{2h)(x)   |   Y(x) - y_{2h}(x)   |   y_{h}(x)   |   Y(x) - y_{h}(x)   |   [y_{h}(x) - y_{2h}(x)]/3   |\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
for i = 1:length(t)
    Y(i) = 1/(1+t(i)); % exact solution of y'
    fprintf('   %1.1f      %1.6f           %1.6f            %1.6f           %1.6f                  %1.6f\n',...
        [t(i) y_2h(i) Y(i)-y_2h(i) y_h(i) Y(i)-y_h(i) (y_h(i)-y_2h(i))/3]);
end
fprintf("--------------------------------------------------------------------------------------------------------------------\n");

%% MATH7003-00: Assignment #3-(2), 2019310290 Sangman Jung
clear,clc
