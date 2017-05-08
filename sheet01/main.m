%% 2. Explicit integrators

%% 2.a
x0 = [1;0]; % this is a column vector
T_s = 0.5;
k = 0:1:19; % to create an array from 0 to 19
t = k*T_s; % This means that the scalar T_s is multiplied with all entries of k
T = t(end);

%[t,x] = ode45(@(t,x) [x(2); -0.2*x(2) - x(1)] , t, x0);
ode = @(t,x)[x(2); -0.2*x(2) - x(1)];
[t,x] = ode45(ode, t, x0);

figure;
hold on;
plot(x(:,1),x(:,2),'-o');

%% 2.b

x_euler = sim_euler(T,0.125,ode,x0);

x_rk4 = rk4(T, 0.5, ode, x0);

plot(x_euler(:,1),x_euler(:,2), '-ro');
plot(x_rk4(:,1),x_rk4(:,2), '-go');
hold off;

%% 2.c

clear all;

import casadi.*;

x = SX.sym('x',1,2);

f = Function('f', {x}, {[x(2), -0.2*x(2) - x(1)]});

x0 = SX([1,2]);

T = SX(19);
h = SX(0.5);

rk4 = Function('rk4', {T, h, x, x0},...
               {git p},...
               {'T', 'h', 'x', 'x0'},...
               {'x_sim'});