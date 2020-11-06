%% MATH7003-00: Assignment #5-(2-2), 2019310290 Sangman Jung.
clear,clc

%%% fundamental setting
h = 0.5; % step size
lambda = [-1 -10 -50]; % parameters
% given the differential equation
f = @(x,y,lambda) lambda*y + (1-lambda)*cos(x) - (1+lambda)*sin(x);
Y = @(x) sin(x)+cos(x); % exact solution of y'

%%% trapezoidal method iteration
for lam_iter = 1:3 % lambda = -1, -10, -50
    y0 = 1; % initial condition
    x0 = 0; % initial x
    x = x0:h:10; % variable x
    num_x = (10-x0)/h;
    y = zeros(1,num_x);
    y(1) = y0;
    for j = 1:num_x % Trapezoidal method loop
        y(j+1) = ((1+(h/2)*lambda(lam_iter))*y(j)+h/2*((1-lambda(lam_iter))*cos(x(j))-...
            (1+lambda(lam_iter))*sin(x(j))+(1-lambda(lam_iter))*cos(x(j+1))-...
            (1+lambda(lam_iter))*sin(x(j+1))))/(1-(h/2)*lambda(lam_iter)); % explicit version
    end
    fval(lam_iter) = {y'}; % save the values of y for each lambda
    xval(lam_iter) = {x'}; % save the values of x for each lambda
end

%%% print out setting
% load the value of x and y for each lambda
x_lam1 = cell2mat(xval(1)); % x for lambda = -1
x_lam2 = cell2mat(xval(2)); % x for lambda = -10
x_lam3 = cell2mat(xval(3)); % x for lambda = -50
y_lam1 = cell2mat(fval(1)); % numerical y for lambda = -1
y_lam2 = cell2mat(fval(2)); % numerical y for lambda = -10
y_lam3 = cell2mat(fval(3)); % numerical y for lambda = -50

%%% table view setting
table = 2:2:10;
for k = table
    index_lam1(k) = find(x_lam1 == k);
    index_lam2(k) = find(x_lam2 == k);
    index_lam3(k) = find(x_lam3 == k);
end
index_lam1 = nonzeros(index_lam1)';
index_lam2 = nonzeros(index_lam2)';
index_lam3 = nonzeros(index_lam3)';

%%% output
Error_lam1 = Y(x_lam1(index_lam1))-y_lam1(index_lam1);
Error_lam2 = Y(x_lam2(index_lam2))-y_lam2(index_lam2);
Error_lam3 = Y(x_lam3(index_lam3))-y_lam3(index_lam3);

%%% print the Table 6.18
fprintf("Table 6.18  Example of trapezoidal rule: h = 0.5\n");
fprintf("--------------------------------------------------------------------------------------\n");
fprintf("|   x   |   Error: lambda = -1   |   Error: lambda = -10   |   Error: lambda = -50   |\n");
fprintf("--------------------------------------------------------------------------------------\n");
for k = 1:length(table)
    fprintf('   %2s           % 1.2e                 % 1.2e               % 1.2e\n',...
        num2str(table(k)), Error_lam1(k), Error_lam2(k), Error_lam3(k));
end
fprintf("--------------------------------------------------------------------------------------\n");
