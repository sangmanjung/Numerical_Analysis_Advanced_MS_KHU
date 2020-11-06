%% MATH7003-00: Assignment #8, 2019310290 Sangman Jung.
clear,clc
% the problem for a linear equation Ax=b
m_iter = 7; % set the iteration numbers
A = [10 3 1; 2 -10 3; 1 3 10]; % matrix A
b = [14 -5 14]'; % vector b
true_x = [1 1 1]'; % solution of Ax=b
x(:,1) = [0 0 0]'; % initial guess for Gauss-Jacobi (GB)
x_sei(:,1) = x(:,1); % initial guess for Gauss-Seidal (GS)
e(:,1) = true_x - x(:,1); % error for (GB)
e_sei(:,1) = true_x - x_sei(:,1); % error for (GS)
norm_e(:,1) = norm(e(:,1),inf); % infinite norm of error 'e' (GB)
norm_e_sei(:,1) = norm(e_sei(:,1),inf); % infinite norm of error 'e' (GS)
Ratio = zeros(m_iter-1,1); % ratio == ||e^{(m)}||_{inf} / ||e^{(m-1)}||_{inf}
Ratio_sei = Ratio;

% Problem 1
fprintf('Table 8.3   Numerical results for the Gauss-Jacobi method\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
fprintf('\tm\t\t{x_{1}}^{(m)}\t\t{x_{2}}^{(m)}\t\t{x_{3}}^{(m)}\t\t||e^{(m)}||_{inf}\t\tRatio\t\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
fprintf('\t%d\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.2f\t\n',0,x(1,1),x(2,1),x(3,1),norm_e(:,1),'');
for m = 1:m_iter-1 % m iteration
    % search (i,j)-entry for the matrix A
    for i = 1:size(A,1) % position i
        % allocations
        sum_ax = 0;
        seisum1 = 0;
        seisum2 = 0;
        for j = 1:size(A,1) % position j
            if i ~= j
                sum_ax = sum_ax + A(i,j)*x(j,m); % inner product for the rows of A and the vector x
            end
        end
        x(i,m+1) = (b(i)-sum_ax)/A(i,i); % compute 'x' using Gauss-Jacobi method
        
        if i-1 == 0
            seisum1 = 0; % consideration for the end is zero
        else
            for j = 1:i-1
                seisum1 = seisum1 + A(i,j)*x_sei(j,m+1); % the first summation in GS
            end
        end
        for j = i+1:size(A,1)
            seisum2 = seisum2 + A(i,j)*x_sei(j,m); % the second summation in GS
        end
        x_sei(i,m+1) = (b(i)-seisum1-seisum2)/A(i,i); % compute 'x' using Gauss-Seidel method
    end
    % error update
    e(:,m+1) = true_x - x(:,m+1);
    e_sei(:,m+1) = true_x - x_sei(:,m+1);
    % infinite norm of error update
    norm_e(:,m+1) = norm(e(:,m+1),inf);
    norm_e_sei(:,m+1) = norm(e_sei(:,m+1),inf);
    % ratio update
    Ratio(m) = norm_e(:,m+1)/norm_e(:,m);
    Ratio_sei(m) = norm_e_sei(:,m+1)/norm_e_sei(:,m);
    fprintf('\t%d\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.2g\t\n',...
        m,x(1,m+1),x(2,m+1),x(3,m+1),norm_e(:,m+1),Ratio(m)); % print the values
end
fprintf('--------------------------------------------------------------------------------------------------------\n');

% Problem 2
fprintf('\n\nTable 8.4   Numerical results for the Gauss-Seidel method\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
fprintf('\tm\t\t{x_{1}}^{(m)}\t\t{x_{2}}^{(m)}\t\t{x_{3}}^{(m)}\t\t||e^{(m)}||_{inf}\t\tRatio\t\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');
fprintf('\t%d\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.2f\t\n',0,x(1,1),x(2,1),x(3,1),norm_e(:,1),'');
for m = 1:m_iter-1
    fprintf('\t%d\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.6f\t\t\t%1.2e\t\t\t%1.2g\t\n',...
        m,x_sei(1,m+1),x_sei(2,m+1),x_sei(3,m+1),norm_e_sei(:,m+1),Ratio_sei(m)); % print the values
end
fprintf('--------------------------------------------------------------------------------------------------------\n');
