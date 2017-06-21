clear all;

import casadi.*;

N = 50;
nx = 2;
nu = 1;
x0 = [pi;0];
rho = 0.01;

x = MX.sym('x', nx); % x = [φ, ω]
u = MX.sym('u', nu); % u = [τ] (torque)

l = (x'*x + 2*u^2 - rho*(log(-u+1)+log(u+1)));
L = Function('L', {x,u}, {l});

ode = [x(2); sin(x(1)) + u];
ODE = Function('ode', {x, u}, {ode});

%%

% w is a combined vector of [x_0(1), x_0(2), u_0,
% x_2(1), ..., u_{N-1},x_N(1), x_N(2)
nw = N*nx + (N-1)*nu;
h = 0.1;

w = MX.sym('w', nw);
r = 0;

G = [];

% fixing x_0
lbx = [pi; 0; -inf];
ubx = [pi; 0; inf];

% for separating x and u later
x_mask = false(nw,1);
u_mask = false(nw,1);

for i=1:N-1
    % helping variable to make formating easier
    j = 3*(i-1) + 1;
    k = 2*(i-1) + 1;
    
    xi = w(j:j+1);
    ui = w(j+2);
    xnext = w(j+3:j+4);
    
    x_mask(j) = true;
    x_mask(j+1) = true;
    u_mask(j+2) = true;

    r = r + L(xi, ui);
    
    % RK4
    k_1 = ODE(xi, ui);
    k_2 = ODE(xi + k_1*h/2, ui);
    k_3 = ODE(xi + k_2*h/2, ui);
    k_4 = ODE(xi + k_3*h, ui);
         
    x_pred = xi + h * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;
    
    % define constraint
    G = [G; x_pred - xnext];
    
    % define upper and lower bounds;
    lbxi = [-pi;-inf;-inf];
    ubxi = [pi;inf;inf];
    
    lbx = [lbx; lbxi];
    ubx = [ubx; ubxi];
end
% remove the last u-bound
lbx = lbx(1: end - 1);
ubx = ubx(1: end - 1);

% add E
r = r + 10* norm(w(end-1:end))^2;

% the last two w coordinates are x as well
x_mask(end-1:end) = [true;true];

%%
j = jacobian(r,w);
J = Function('J',{w},{j});

h = hessian(r,w);
H = Function('H',{w},{h});

hessian_LaGrange = full(H(0.1*ones(nw,1)));
%spy(hessian_LaGrange);
%%
nlp = struct('x', w, 'f', r, 'g', G);

solver = nlpsol('solver', 'ipopt', nlp);

res = solver('x0' , [0.1*ones(nw, 1)],...       % solution guess
             'lbx', lbx,...  % lower bound on x
             'ubx', ubx,...  % upper bound on x
             'lbg', [zeros(size(G))],...
             'ubg', [zeros(size(G))]);

%% reshape the w-vector

result = full(res.x);
x_res_sim = zeros(N,nx);
dummy = result(x_mask);
% gets every second entry
x_res_sim(:,1) = dummy(1:2:length(dummy));
x_res_sim(:,2) = dummy(2:2:length(dummy));
%x_res_sim = reshape(result(x_mask),[],2);
u_res_sim = result(u_mask);

%% plot the stuff!

control_fig = figure;
hold on;
plot(1:N-1,u_res_sim,'r-');

title('Control');
legend('u');
hold off;

trajectory_fig = figure;
hold on;
plot(1:N,x_res_sim(:,1));
plot(1:N,x_res_sim(:,2));
title('Trajectory');
legend('φ', 'ω');

hold off;