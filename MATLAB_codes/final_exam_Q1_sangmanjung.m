%% final exam: problem #Q1. 2019310290 Sangman Jung

clear,clc
N = 16; h = 1/N; delta = h; % step N and step size h and delta
x = 0:h:1; y = 0:delta:1; % spatial variables
nx = length(x); ny = length(y); % the size of x and y
ksi = 1-2*((sin(pi/(2*N)))^2);
omega = 2/(1+sqrt(1-(ksi)^2)); % optimal acceleration parameter
epsilon = 0.001; % stop criterion
f = @(x,y) (x^2)*y^2; % exact solution of our problem
g = @(x,y) 2*(x^2+y^2); % choose the function g

for w_iter = 1:3  % w_iter = 1 is G-J, w_iter = 2 is G-S, w_iter = 3 is SOR
    u_ini = zeros(nx,ny); % initial guess
    u_exact = zeros(nx,ny); % pre-allocation of exact solution
    G = zeros(nx,ny); % pre-allocation of the function g
    u_new = zeros(nx,ny); % pre-allocation of numerical solution
    m_iter = 1; % initialize the iteration number
    
    for k = 1:ny
        for j = 1:nx
            u_ini(j,k)=(1-x(j))*f(0,y(k))+x(j)*f(1,y(k))+(1-y(k))*f(x(j),0)+...
                y(k)*f(x(j),1)-((1-y(k))*(1-x(j))*f(0,0)+(1-y(k))*x(j)*f(1,0)+...
                y(k)*(1-x(j))*f(0,1)+x(j)*y(k)*f(1,1)); % compute the initial guess
            u_exact(j,k)=f(x(j),y(k)); % allocate the exact solution
            G(j,k)=g(x(j),y(k)); % allocate the function g
        end
    end
    u_old = u_ini; % update the m-iteration
    
    % boundary values of the problem
    u_new(1,:) = u_exact(1,:); % left vertical
    u_new(:,1) = u_exact(:,1); % bottom
    u_new(end,:) = u_exact(end,:); % right vertical
    u_new(:,end) = u_exact(:,end); % top
    
    error = 1; % initialize the error
    while 1 % Main Loop
        for k = 2:ny-1
            for j = 2:nx-1
                if w_iter == 1
                    % Gauss-Jacobi method (G-J)
                    u_new(j,k) = (u_old(j+1,k)+u_old(j,k+1)+u_old(j-1,k)+u_old(j,k-1))/4-(h^2)*G(j,k)/4;
                elseif w_iter > 1
                    % Gauss-Seidel method (G-S)
                    u_new(j,k) = (u_old(j+1,k)+u_old(j,k+1)+u_new(j-1,k)+u_new(j,k-1))/4-(h^2)*G(j,k)/4;
                    if w_iter > 2
                        % Successive OverRelaxation (SOR) method
                        u_new(j,k) = omega*u_new(j,k) + (1-omega)*u_old(j,k);
                    end
                end
            end
        end
        if m_iter > 1 % compute the error (8.7.5) in Chapter 8.7
            c = max(max(abs(u_new-u_old)))/max(max(abs(u_old-u_2old)));
            error = c/(1-c)*(max(max(abs(u_new-u_old))));
            actual_error(m_iter,w_iter) = max(max(abs(u_exact-u_new)));
            iter_error(m_iter,w_iter) = max(max(abs(u_new-u_old)));
            est_c(m_iter,w_iter) = c;
            est_error(m_iter,w_iter) = error; 
        end
        if error <= epsilon % error criterion
            break;
        end
        m_iter = m_iter + 1; % update the error after passing the criterion
        u_2old = u_old; % uptate from m-1 to m
        u_old = u_new; % update from m to m+1
    end
    uval(w_iter) = {u_new}; % save the numerical solution
end

for i = 1:3
    itm(i) = size(find(est_c(:,i) ~=0),1);
end

fprintf('\n Q1.(a): The results of the Table for Gauss-Jacobi method\n')
fprintf('       (Let N = 16, epsilon = 0.001, max_iter = %d) \n',itm(1))
fprintf('-----------------------------------------------------------\n')
fprintf('      m     Iter_Error   Est_Error     Error       c       \n')
fprintf('-----------------------------------------------------------\n')
for t = 20:25
    fprintf('    %3.0f      %1.2e    %1.2e    %1.2e    %3.3f        \n',...
        [t  iter_error(t,1) est_error(t,1) actual_error(t,1) est_c(t,1)]);
end
fprintf('-----------------------------------------------------------\n')

fprintf('\n       Q1.(b): The results of the Table for Gauss-Seidel method     \n')
fprintf('             (Let N = 16, epsilon = 0.001, max_iter = %d)      \n',itm(2))
fprintf('----------------------------------------------------------------------\n')
fprintf('      m     Iter_Error   Est_Error     Error   Est_Ratio   Ratio \n')
fprintf('----------------------------------------------------------------------\n')
for t = 20:25
    fprintf('    %3.0f      %1.2e    %1.2e    %1.2e    %3.3f     %3.3f         \n',...
        [t  iter_error(t,2) est_error(t,2) actual_error(t,2) est_c(t,2) 1-(pi*h)^2]);
end
fprintf('----------------------------------------------------------------------\n')

fprintf('\n          Q1.(c): The results of the Table for SOR method              \n')
fprintf('            (Let N = 16, epsilon = 0.001, max_iter = %d)          \n',itm(3))
fprintf('----------------------------------------------------------------------\n')
fprintf('      m     Iter_Error   Est_Error     Error   Est_Ratio   Ratio \n')
fprintf('----------------------------------------------------------------------\n')
for t = 14:19
    fprintf('    %3.0f      %1.2e    %1.2e    %1.2e    %3.3f     %3.3f         \n',...
        [t  iter_error(t,3) est_error(t,3) actual_error(t,3) est_c(t,3) omega-1]);
end
fprintf('----------------------------------------------------------------------\n')