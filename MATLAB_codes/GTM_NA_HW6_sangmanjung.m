%% MATH7003-00: Assignment #6-1, 2019310290 Sangman Jung.
%%% Method of Lines: Euler's method %%%
clear,clc
% exact solution U and the function G
U = @(x,t) exp(-0.1*t).*sin(pi*x);
G = @(x,t) U(x,t)*(pi^2-0.1);
% set the parameters
m = [4 8 16]; % the length of spatial variable x
delta = 1./m; % space step delta
% method of lines loop
for m_iter = 1:length(m) % m = 4, 8, 16
    x = (0:m(m_iter))*delta(m_iter); % calculate the variable x
    h = (delta(m_iter)^2)/2; % step size depend on delta
    t = 0:h:5; % time variable t
    A = zeros(m(m_iter)-1,m(m_iter)-1); % allocate the matrix 'LAMBDA'
    A = (1/delta(m_iter)^2)*(diag(diag(-2*ones(size(A)))) +...
        diag(ones(size(A,1)-1,1),-1) +...
        diag(ones(size(A,1)-1,1),1)); % the matrix LAMBDA
    d0 = U(x(1),t); d1 = U(x(end),t); % boundary value of U
    u0 = U(x(2:end-1),0)'; % initial value of U
    g = zeros(length(x)-2,length(t)); % allocate the function g
    g(1,:) = (d0/delta(m_iter)^2)+G(x(2),t);
    g(end,:) = (d1/delta(m_iter)^2)+G(x(end-1),t);
    g(2:end-1,:) = G(x(3:end-2)',t); % the function g
    % Euler's method
    V = zeros(length(x)-2,1); % allocate the value of V
    V(:,1) = u0; % initial condition for the method
    for n = 1:length(t)-1
        V(:,n+1) = V(:,n)+h*(A*V(:,n)+g(:,n)); % Euler's method
    end
    u = zeros(length(x),length(t)); % allocate the value of U
    u(1,:) = d0; u(end,:) = d1; % boundary condition
    u(2:end-1,:) = V; % save the value of V
    % recode the values
    xval(m_iter) = {x'};
    tval(m_iter) = {t};
    Uval(m_iter) = {u};
    Vval(m_iter) = {V};
end
% setting the values for table 6.20
x4 = cell2mat(xval(1));
t4 = cell2mat(tval(1));
u4 = cell2mat(Uval(1));
x8 = cell2mat(xval(2));
t8 = cell2mat(tval(2));
u8 = cell2mat(Uval(2));
x16 = cell2mat(xval(3));
t16 = cell2mat(tval(3));
u16 = cell2mat(Uval(3));
Table = 1:5;
for i = 1:5
    ind_t4(i) = find(t4 == i);
    ind_t8(i) = find(t8 == i);
    ind_t16(i) = find(t16 == i);
end
Error4 = max(abs(U(x4,t4(ind_t4))-u4(:,ind_t4)));
Error8 = max(abs(U(x8,t8(ind_t8))-u8(:,ind_t8)));
Error16 = max(abs(U(x16,t16(ind_t16))-u16(:,ind_t16)));
% Table 6.20
fprintf("Table 6.20 The method of lines: Euler's method\n")
fprintf('-----------------------------------------------------------\n')
fprintf('          Error                Error               Error \n')
fprintf('  t       m = 4     Ratio      m = 8     Ratio     m = 16\n')
fprintf('-----------------------------------------------------------\n')
for i = 1:5
    fprintf('% 1.1f   % 1.2e   % 1.2f   % 1.2e    % 1.2f    % 1.2e\n',...
        [Table(i) Error4(i) Error4(i)/Error8(i) Error8(i) Error8(i)/Error16(i) Error16(i)]);
end
fprintf('-----------------------------------------------------------\n\n\n')

%% MATH7003-00: Assignment #6-2, 2019310290 Sangman Jung.
%%% Method of Lines: Backward Euler method %%%
clear
% exact solution U and the function G
U = @(x,t) exp(-0.1*t).*sin(pi*x);
G = @(x,t) U(x,t)*(pi^2-0.1);
% set the parameters
m = [4 8 16]; % the length of spatial variable x
delta = 1./m; % space step delta
% method of lines loop
for m_iter = 1:length(m) % m = 4, 8, 16
    x = (0:m(m_iter))*delta(m_iter); % calculate the variable x
    h = 0.1; % step size not depend on delta
    t = 0:h:5; % time variable t
    A = zeros(m(m_iter)-1,m(m_iter)-1); % allocate the matrix 'LAMBDA'
    A = (1/delta(m_iter)^2)*(diag(diag(-2*ones(size(A)))) +...
        diag(ones(size(A,1)-1,1),-1) +...
        diag(ones(size(A,1)-1,1),1)); % the matrix LAMBDA
    d0 = U(x(1),t); d1 = U(x(end),t); % boundary value of U
    u0 = U(x(2:end-1),0)'; % initial value of U
    g = zeros(length(x)-2,length(t)); % allocate the function g
    g(1,:) = (d0/delta(m_iter)^2)+G(x(2),t);
    g(end,:) = (d1/delta(m_iter)^2)+G(x(end-1),t);
    g(2:end-1,:) = G(x(3:end-2)',t); % the function g
    % backward Euler method
    I = diag(diag(ones(size(A)))); % identity matrix
    V = zeros(length(x)-2,1); % allocate the value of V
    V(:,1) = u0; % initial condition for the method
    for n = 1:length(t)-1
        V(:,n+1) = inv(I-h*A)*(V(:,n)+h*g(:,n+1)); % backward Euler method
    end
    u = zeros(length(x),length(t)); % allocate the value of U
    u(1,:) = d0; u(end,:) = d1; % boundary condition
    u(2:end-1,:) = V; % save the value of V
    % recode the values
    xval(m_iter) = {x'};
    tval(m_iter) = {t};
    Uval(m_iter) = {u};
    Vval(m_iter) = {V};
end
% setting the values for table 6.20
x4 = cell2mat(xval(1));
t4 = cell2mat(tval(1));
u4 = cell2mat(Uval(1));
x8 = cell2mat(xval(2));
t8 = cell2mat(tval(2));
u8 = cell2mat(Uval(2));
x16 = cell2mat(xval(3));
t16 = cell2mat(tval(3));
u16 = cell2mat(Uval(3));
Table = 1:5;
for i = 1:5
    ind_t4(i) = find(t4 == i);
    ind_t8(i) = find(t8 == i);
    ind_t16(i) = find(t16 == i);
end
Error4 = max(abs(U(x4,t4(ind_t4))-u4(:,ind_t4)));
Error8 = max(abs(U(x8,t8(ind_t8))-u8(:,ind_t8)));
Error16 = max(abs(U(x16,t16(ind_t16))-u16(:,ind_t16)));
% Table 6.21
fprintf("Table 6.21 The method of lines: backward Euler method\n")
fprintf('-------------------------------------------------------\n')
fprintf('             Error            Error           Error \n')
fprintf('    t        m = 4            m = 8           m = 16\n')
fprintf('-------------------------------------------------------\n')
for i = 1:5
    fprintf('  % 1.1f    % 1.2e        % 1.2e        % 1.2e\n',...
        [Table(i) Error4(i) Error8(i) Error16(i)]);
end
fprintf('-------------------------------------------------------\n')
