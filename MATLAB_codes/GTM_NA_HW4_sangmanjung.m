%% Homework # 4, 2019310290 Sangman Jung
clear,clc
% define ODE and it's exact solution
f = @(x,y) (1/(1+x^2)) - 2*y^2; % ODE in the example
Y = @(x) x/(1+x^2); % exact solution of y'
% initial values
y0 = 0; % y(0) = 0
% interval of x
x0 = 0; x_end = 10; % x in [0,10]
% step size h
h_min = 0.001; h_max = 1; h = 0.1;
% error tolerence
epsilon = 0.0005;
% run Detrap
DetrapHW4(Y,f,x0,y0,x_end,epsilon,h,h_min,h_max)