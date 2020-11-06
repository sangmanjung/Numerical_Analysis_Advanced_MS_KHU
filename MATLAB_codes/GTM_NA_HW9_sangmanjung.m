%% MATH7003-00: Assignment #9, 2019310290 Sangman Jung.
clear,clc

N = [8 16 32]; % step N
epsilon = [0.01 0.001 0.0001]; % error criterion
f = @(x,y) exp(pi*x).*cos(pi*y); % exact solution of our problem

ksi = 1-2.*((sin(pi./(2*N))).^2);
omega = 2./(1+sqrt(1-(ksi).^2)); % optimal acceleration parameter

uval{length(N),length(epsilon),2} = []; % pre-allocation of the numerical solution

iteration = zeros(length(N),length(epsilon)); % pre-allocation of G-S iteration
iteration_sor = zeros(length(N),length(epsilon)); % pre-allocation of SOR iteration 

% Full Loop
for w_iter = 1:2 % w_iter = 1 is G-S, w_iter = 2 is SOR
    for N_iter = 1:length(N) % iterate per step N
        h = 1/N(N_iter); % step size
        x = 0:h:1; y = x; % spatial variables (grid)
        nx = length(x); ny = length(y); % the size of the grid
        u_ini = zeros(nx,ny); % initial guess
        u_exact = zeros(nx,ny); % pre-allocation of exact solution
        u_new = zeros(nx,ny); % pre-allocation of numerical solution
        for eps_iter = 1:length(epsilon) % iteration of the criterion
            m_iter = 1; % initialize the iteration number
            for k = 1:ny
                for j = 1:nx
                    u_ini(j,k)=(1-x(j))*f(0,y(k))+x(j)*f(1,y(k))+(1-y(k))*f(x(j),0)+...
                        y(k)*f(x(j),1)-((1-y(k))*(1-x(j))*f(0,0)+(1-y(k))*x(j)*f(1,0)+...
                        y(k)*(1-x(j))*f(0,1)+x(j)*y(k)*f(1,1)); % compute the initial guess
                    u_exact(j,k)=f(x(j),y(k)); % allocate the exact solution
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
                        % Gauss-Seidel (G-S) method
                        u_new(j,k) = (u_old(j+1,k)+u_old(j,k+1)+u_new(j-1,k)+u_new(j,k-1))/4;
                        if w_iter > 1
                            % Successive OverRelaxation (SOR) method
                            u_new(j,k) = omega(N_iter)*u_new(j,k) + (1-omega(N_iter))*u_old(j,k);
                        end
                    end
                end
                if m_iter > 1 % compute the error (8.7.5) in Chapter 8.7
                    c = max(max(abs(u_new-u_old)))/max(max(abs(u_old-u_2old)));
                    error = c/(1-c)*(max(max(abs(u_new-u_old))));
                end
                if error <= epsilon(eps_iter) % error criterion
                    break;
                end
                m_iter = m_iter + 1; % update the error after passing the criterion
                u_2old = u_old; % uptate from m-1 to m
                u_old = u_new; % update from m to m+1
            end
            if w_iter == 1 % G-S iteration
                iteration(N_iter,eps_iter) = m_iter;
            else % SOR iteration
                iteration_sor(N_iter,eps_iter) = m_iter;
            end
            uval(eps_iter,N_iter,w_iter) = {u_new}; % save the numerical solution
        end
    end
end

% summarize our results in order to obtain the table
N_table = [N(1) N N]';
eps_table = [epsilon(1) repmat(epsilon(2),1,3) repmat(epsilon(3),1,3)]';
iteration = reshape(iteration,[9 1]);
iteration = iteration([1 4:9]);
iteration_sor = reshape(iteration_sor,[9 1]);
iteration_sor = iteration_sor([1 4:9]);

% print Table 8.9
fprintf('\nTable 8.9  Number of iterates necessary to solve (8.8.5)\n')
fprintf('---------------------------------------------------------\n')
fprintf('     N        epsilon        Gauss-Seidel        SOR     \n')
fprintf('---------------------------------------------------------\n')
for t = 1:length(N_table)
    fprintf('   %3.0f        %6.4g             %3.0f             %3.0f\n',...
        [N_table(t) eps_table(t) iteration(t) iteration_sor(t)]);
end
fprintf('---------------------------------------------------------\n')
