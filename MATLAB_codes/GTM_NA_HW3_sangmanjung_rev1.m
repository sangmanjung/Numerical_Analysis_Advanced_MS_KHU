%% MATH7003-00: Assignment #3-(1), 2019310290 Sangman Jung _ revision 1.
clear,clc

% Initial parameters, initial conditions
h = [0.25 0.5]; % step size
y(1) = 1; % initial value of y
Y = @(x) 1./(1+x)'; % exact solution of y'

% Trapezoidal method iteration to solve y'=-y^2
% Using predictor-corrector method, we can obtain the scheme as follows.
for h_iter = 1:2 % the loop of h and 2h
    x = 0:h(h_iter):5; % apply a different step size
    y = ones(size(x)); % reset the computed results of h=0.25 before
    y(2) = y(1)-h(h_iter)*y(1)^2; % define the term for midpoint, using Euler's method
    y_p(1)=y(1)-2*h(h_iter)*y(2)^2; % initial guess of y_p (midpoint method)
    for j =1:length(x)-1 % Trapezoidal method loop
        y_p(j+1)=y(j)-h(h_iter)*(y(j)^2+y_p(j)^2)/2; % predictor
        y(j+1)=y(j)-h(h_iter)*(y(j)^2+y_p(j+1)^2)/2; % corrector (In this code, numerical sol. of y')
    end
    fval(h_iter) = {y'}; % save the values of y have different size for y_h and y_2h
    xval(h_iter) = {x'}; % save the values of x have different size for y_h and y_2h
end

% Print out setting
x_h = cell2mat(xval(1)); % for the case of h = 0.25
x_2h = cell2mat(xval(2)); % for the case of h = 0.5
y_h = cell2mat(fval(1)); % numerical y for h = 0.25
y_2h = cell2mat(fval(2)); % numerical y for h = 0.5
y_ht = y_h(5:4:end);
y_2ht = y_2h(3:2:end);
t = 1:5;

% Print the Table 6.8
fprintf("Table 6.8  Trapezoidal method and Richardson error estimation\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
fprintf("|   x   |   y_{2h)(x)   |   Y(x) - y_{2h}(x)   |   y_{h}(x)   |   Y(x) - y_{h}(x)   |   [y_{h}(x) - y_{2h}(x)]/3   |\n");
fprintf("--------------------------------------------------------------------------------------------------------------------\n");
for i = 1:length(t)
    fprintf('   %1.1f      %1.6f           %1.6f            %1.6f           %1.6f                  %1.6f\n',...
        abs([t(i) y_2ht(i) Y(t(i))-y_2ht(i) y_ht(i) Y(t(i))-y_ht(i) (y_ht(i)-y_2ht(i))/3]));
end
fprintf("--------------------------------------------------------------------------------------------------------------------\n");

%% Graphs for the results
subplot(1,2,1);
plot(x_h,y_h,'-or',x_2h,y_2h,'-vb',x_h,1./(1+(0:0.25:5)),'-sk');
xlim([0 5]);
xlabel('x-axis');
ylabel('solutions for  {dy/dx = -y^2}');
legend('y_h(x)','y_{2h}(x)','Y(x)');
title('A comparison of each solution');
tt = get(gca,'title');
tt.FontWeight = 'bold';

for h_iter = 1:2 % the loop of h and 2h
    x = 0:h(h_iter):31; % apply a different step size
    y = ones(size(x)); % reset the computed results of h=0.25 before
    y(2) = y(1)-h(h_iter)*y(1)^2; % define the term for midpoint, using Euler's method
    y_p(1)=y(1)-2*h(h_iter)*y(2)^2; % initial guess of y_p (midpoint method)
    for j =1:length(x)-1 % Trapezoidal method loop
        y_p(j+1)=y(j)-h(h_iter)*(y(j)^2+y_p(j)^2)/2; % predictor
        y(j+1)=y(j)-h(h_iter)*(y(j)^2+y_p(j+1)^2)/2; % corrector (In this code, numerical sol. of y')
    end
    fval(h_iter) = {y'}; % save the values of y have different size for y_h and y_2h
    xval(h_iter) = {x'}; % save the values of x have different size for y_h and y_2h
end

% Print out setting
x_h = cell2mat(xval(1)); % for the case of h = 0.25
x_2h = cell2mat(xval(2)); % for the case of h = 0.5
y_h = cell2mat(fval(1)); % numerical y for h = 0.25
y_2h = cell2mat(fval(2)); % numerical y for h = 0.5
y_ht = y_h(5:4:end);
y_2ht = y_2h(3:2:end);
t = 1:31;

% Graphs
subplot(1,2,2);
plot(t,Y(t)-y_ht,'-ro',t,(y_ht-y_2ht)/3,'-bo',t,(Y(t)-y_ht)-(y_ht-y_2ht)/3,'-ko');
xlim([1 31]);
ylim([-0.00055 0.00095])
xticks(1:5:31);
xticklabels({'1','5','10','15','20','25','30'});
xlabel('x-axis')
ylabel('Error estimations')
legend('True error','Richardson error','True - Richardson');
title('A comparison of errors');
tt = get(gca,'title');
tt.FontWeight = 'bold';