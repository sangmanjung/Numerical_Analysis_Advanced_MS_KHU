%% final exam: problem #Q2. 2019310290 Sangman Jung

clear,clc
% Example for chapter 8.9
n = 40;
A = [5 4 3 2 1; 4 5 4 3 2; 3 4 5 4 3; 2 3 4 5 4; 1 2 3 4 5];
b=[7.9380 12.9763 17.3057 19.4332 18.4196]';
true_x = [-0.3227 0.3544 1.1010 1.5705 1.6897]';

% Conjugate Gradient Method
[X,N,r_norm_inf,x_norm_inf,x_norm_A, error_bound] = CG(A,b,true_x,n);

% Table 8.10
fprintf('\n Table 8.10  Example of the conjugate gradient method \n')
fprintf('----------------------------------------------------------------------------------------------------\n')
fprintf('      k     || r_{k} ||_{inf}   || x - x_{k} ||_{inf}   || x - x_{k} ||_{A}      Bound (8.9.20)     \n')
fprintf('----------------------------------------------------------------------------------------------------\n')
for t = 1:5
    fprintf('      %d          %1.2e              %1.2e               %1.2e                %2.2f        \n',...
        [N(t) r_norm_inf(t) x_norm_inf(t) x_norm_A(t) error_bound(t)]);
end
fprintf('----------------------------------------------------------------------------------------------------\n')

% Additional Table for the true solution and numerical solution
fprintf('\nTwo solutions in Table 8.10 \n')
fprintf('   (max_iteration = %d)     \n',n)
fprintf('----------------------------\n')
fprintf('     x^{*}    numerical x   \n')
fprintf('----------------------------\n')
for t = 1:length(true_x)
    fprintf('    % 3.4f    % 1.4f   \n',[true_x(t) X(t)])
end
fprintf('----------------------------\n')