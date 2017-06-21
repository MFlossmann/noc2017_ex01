clear all;

import casadi.*

%% Task 1 - Sequential approach

N = 50;
nx = 2;
nu = 1;
x0 = [pi;0];

x = MX.sym('x', nx); % x = [φ, ω]
u = MX.sym('u', nu); % u = [τ] (torque)

l = (x'*x + 2*u^2);
L = Function('L', {x,u}, {l});

ode = [x(2), sin(x(1)) + u];
ODE = Function('ode', {x, u}, {ode});

%%

% RK4 simulation
x_bar = MX(N, nx);
U = MX.sym('U',N-1, nu);

x_bar(1,:) = x0;
h = 0.1;
phi = 0;
for i=1:N-1    
    k_1 = ODE(x_bar(i,:), U(i));
    k_2 = ODE(x_bar(i,:) + k_1*h/2, U(i));
    k_3 = ODE(x_bar(i,:) + k_2*h/2, U(i));
    k_4 = ODE(x_bar(i,:) + k_3*h, U(i));
         
    x_bar(i+1,:) = x_bar(i,:) + h * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;
    
    phi = phi + L(x_bar(i,:), U(i));
end

X = Function('X', {U}, {x_bar});

% +E(x)
phi = phi + 10* norm(x_bar(end))^2;

Phi = Function('Phi', {U}, {phi});

h = hessian(Phi(U),U);
H = Function('H', {U}, {h});

U_init = 0.1*ones(N-1,1);
spy(full(H(U_init)));
%% Solver

nlp = struct('x', U, 'f', phi);

solver = nlpsol('solver', 'ipopt', nlp);

res = solver('x0' , U_init,...       % solution guess
             'lbx', -inf*ones(N-1,1),...  % lower bound on x
             'ubx', inf*ones(N-1,1));  % upper bound on x

u_res_seq = full(res.x);

%% Plot the trajectories

x_res_seq = full(X(u_res_seq));
trajectory_fig = figure;
hold on;
plot(1:N, x_res_seq(:,1));
plot(1:N, x_res_seq(:,2));
title('Trajectory');
legend('φ_{seq}', 'ω_{seq}');
hold off;

control_fig = figure;
hold on;
plot(1:N-1, u_res_seq);
title('Control');
legend('τ_seq');
hold off;

%%
clearvars -except trajectory_fig control_fig u_res_seq x_res_seq