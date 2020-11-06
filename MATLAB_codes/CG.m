function [X,N,r_norm_inf,x_norm_inf,x_norm_A, error_bound] = CG(A,b,x,n)
% final exam: problem #Q2. 2019310290 Sangman Jung
% Remark: This algorithm calculates the solution of Ax = b
% using the conjugate gradient method.

% list of iteration
N = 1:n;

% pre-allocation
x0 = zeros(size(x)); r0 = b; p0 = zeros(size(x));
error_bound = zeros(n,1);
r_norm_inf = zeros(n,1);
x_norm_inf = zeros(n,1);
x_norm_A = zeros(n,1);

% compute the eigenvalues of A and c
eigen_value = zeros(size(x));
eigen_value(:) = eig(A);
c = eigen_value(1)/eigen_value(end);

% step 2
x_old = x0; r_old = r0; p_old = p0;

% step 3
for k=1:n-1
    % step 4
    if r_old == 0
        X = x_old;
    end
    % step 5
    if k == 1
        beta = 0;
        r_2old = r_old;
    else
        beta = (r_old'*r_old)/(r_2old'*r_2old);
    end
    p_new = r_old + beta*p_old;
    % step 6
    alp = (r_old'*r_old)/(p_new'*A*p_new);
    x_new = x_old + alp*p_new;
    r_new = b - A*x_new;
    
    % compute the column of Table 8.10
    difference = x - x_new;
    r_norm_inf(k,1) = max(abs(r_new));
    x_norm_inf(k,1) = max(abs(difference));
    x_norm_A(k,1) = sqrt(difference'*A*difference);
    error_bound(k,1) = 2*(((1-sqrt(c))/(1+sqrt(c)))^k)*sqrt(x'*A*x);
    
    % update the results for each k
    x_old = x_new;
    r_old = r_new;
    p_old = p_new;
end

% step 8
if n > 1
    X = x_new;
else
    X = [];
    disp('<Warning> Enter the value n > 1.');
end