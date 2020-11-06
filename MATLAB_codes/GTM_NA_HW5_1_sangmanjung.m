%% MATH7003-00: Assignment #5-(1), 2019310290 Sangman Jung
clear,clc

%%% fundamental setting
h = [0.125 0.0625]; % step size
f = @(x,y) (1./(1+x.^2))-2*y.^2; % differential equation
Y = @(x) x./(1+x.^2); % exact solution
x0 = 0; x_end = 10; y0 = 0; % interval and initial condition
fval = cell(1,length(h)); % save the values of y for y_0125 and y_00625
xval = cell(1,length(h)); % save the values of x for y_0125 and y_00625
num_x = zeros(1,2); % allocation of size of interval x

%%% Adams method iteration: using predictor-corrector method
for h_iter = 1:length(h)
    x = x0:h(h_iter):x_end; % interval x
    num_x(h_iter) = (x_end-x0)/h(h_iter); % interval size for each h
    y = zeros(1,num_x(h_iter)); % allocation of the solution y
    y_p = zeros(1,num_x(h_iter)); % allocation of the predictor
    % The initial values y1,y2,y3 were taken to be the true values to
    % simplify the example. <---- refer to the textbook in our class.
    y(1) = y0; y(2:3) = Y(x(2:3));
    % Adam-Bashforth predictor & Adam-Moulton corrector iteration
    for j = 3:num_x(h_iter)
        % predictor
        y(j+1) = y(j) + (h(h_iter)/12)*(23*f(x(j),y(j))-16*f(x(j-1),y(j-1))+5*f(x(j-2),y(j-2)));
        % corrector
        y(j+1) = y(j) + (h(h_iter)/24)*(9*f(x(j+1),y(j+1))+19*f(x(j),y(j))-5*f(x(j-1),y(j-1)) + f(x(j-2),y(j-2)));
    end
    % save the values
    fval(h_iter) = {y'};
    xval(h_iter) = {x'};
end

%%% print out setting
% numerical x, y for h = 0.125
x_0125 = cell2mat(xval(1));
y_0125 = cell2mat(fval(1));
% numerical x, y for h = 0.0625
x_00625 = cell2mat(xval(2));
y_00625 = cell2mat(fval(2));

% we need to represent the values in Table 6.14.
table = 2:2:10;
index_0125 = zeros(1,length(table)*2);
index_00625 = zeros(1,length(table)*2);
for k = table
    index_0125(k) = find(x_0125 == k);
    index_00625(k) = find(x_00625 == k);
end
index_0125 = nonzeros(index_0125)';
index_00625 = nonzeros(index_00625)';

% the error for each step size and the ratio of the error with two step sizes.
Error_125 = Y(x_0125(index_0125))-y_0125(index_0125);
Error_0625 = Y(x_00625(index_00625))-y_00625(index_00625);
Ratio = Error_125./Error_0625;

% print the Table 6.14
fprintf("Table 6.14  Numerical example of the Adams method\n");
fprintf("--------------------------------------------------------------------------\n");
fprintf("|   x   |   Error for h = 0.125   |   Error for h = 0.0625   |   Ratio   |\n");
fprintf("--------------------------------------------------------------------------\n");
for i = 1:length(table)
    fprintf('   %1.1f            %1.2e                 %1.2e              %2.1f    \n',...
        abs([table(i) Error_125(i) Error_0625(i) Ratio(i)]));
end
fprintf("--------------------------------------------------------------------------\n");
