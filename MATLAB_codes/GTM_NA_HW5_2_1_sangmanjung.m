%% MATH7003-00: Assignment #5-(2-1), 2019310290 Sangman Jung.
clear,clc

%%% fundamental setting
h = [0.5 0.1 0.01]; % step size
lambda = [-1 -10 -50]; % parameters
% given the differential equation
f = @(x,y,lambda) lambda*y + (1-lambda)*cos(x) - (1+lambda)*sin(x);
Y = @(x) sin(x)+cos(x); % exact solution of y'

%%% Euler's method for the ODE :
% y' = lambda*y + (1-lambda)*cos(x) - (1+lambda)*sin(x);
for lam_iter = 1:3 % lambda : -1, -10, -50
    for h_iter = 1:3 % step size : 0.5, 0.1, 0.01
        y = []; % initialize
        y0 = 1; % initial condition
        x0 = 0; % initial x
        x = x0:h(h_iter):10; % variable x
        for j = 1:length(x)-1 % Euler's method loop
            y(j+1) = y0 + h(h_iter)*f(x0,y0,lambda(lam_iter)); % numerical scheme
            y0 = y(j+1); % update y
            x0 = x0 + h(h_iter); % update x per step size
        end
        fval(lam_iter,h_iter) = {y'}; % save the values of y for each lambda and step size
        xval(lam_iter,h_iter) = {x'}; % save the values of for each step size
    end
end

%%% print out setting
% load the value of x and y for each lambda and step size
x_lam1(:,[1 2 3]) = cell2mat(xval([1 2 3]));
x_lam2(:,[1 2 3]) = cell2mat(xval([4 5 6]));
x_lam3(:,[1 2 3]) = cell2mat(xval([7 8 9]));
y_lam1(:,[1 2 3]) = cell2mat(fval([1 2 3]));
y_lam2(:,[1 2 3]) = cell2mat(fval([4 5 6]));
y_lam3(:,[1 2 3]) = cell2mat(fval([7 8 9]));

% table view setting
table = 1:5;
for k = table
    index_h1(k) = find(x_lam1(:,1) == k);
    index_h2(k) = find(x_lam2(:,1) == k);
    index_h3(k) = find(x_lam3(:,1) == k);
end
index_h1 = nonzeros(index_h1)';
index_h2 = nonzeros(index_h2)';
index_h3 = nonzeros(index_h3)';

%%% output
Error_lam1 = Y(x_lam1(index_h1,:))-y_lam1(index_h1,:);
Error_lam2 = Y(x_lam2(index_h2,:))-y_lam2(index_h2,:);
Error_lam3 = Y(x_lam3(index_h3,:))-y_lam3(index_h3,:);

%%% print the Table 6.17
fprintf("Table 6.17  Euler's method for (6.8.51)\n");
fprintf("--------------------------------------------------------------------------------------\n");
fprintf("|   lambda   |   x   |   Error: h = 0.5   |   Error: h = 0.1   |   Error: h = 0.01   |\n");
fprintf("--------------------------------------------------------------------------------------\n");
for k = 1:length(table)
    fprintf('      %d         %d          %+1.2e           %+1.2e             %+1.2e\n',...
        [lambda(1) table(k) Error_lam1(k,1) Error_lam2(k,1) Error_lam3(k,1)]);
end
fprintf("--------------------------------------------------------------------------------------\n");
for k = 1:length(table)
    fprintf('     %d         %d          %+1.2e           %+1.2e             %+1.2e\n',...
        [lambda(2) table(k) Error_lam1(k,2) Error_lam2(k,2) Error_lam3(k,2)]);
end
fprintf("--------------------------------------------------------------------------------------\n");
for k = 1:length(table)
    fprintf('     %d         %d          %+1.2e           %+1.2e             %+1.2e\n',...
        [lambda(3) table(k) Error_lam1(k,3) Error_lam2(k,3) Error_lam3(k,3)]);
end
fprintf("--------------------------------------------------------------------------------------\n");
