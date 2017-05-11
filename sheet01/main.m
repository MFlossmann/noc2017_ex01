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

T = 19;
DT = 0.5;
N = T/DT;

x = MX.sym('x',1,2);

xdot = [x(2); -0.2*x(2) - x(1)]; % symbolic representation of xdot

% Runge-Kutta4

   
   f = Function('f', {x}, {xdot}); % function for applying xdot
   X0 = MX.sym('X0', 2);
   X = MX.sym('X',N,2);
   X(1,:) = X0;
   %%
   for i=2:N
       k1 = f(X(i,:));
       k2 = f(X(i,:) + DT/2 * k1);
       k3 = f(X(i,:) + DT/2 * k2);
       k4 = f(X(i,:) + DT * k3);
       X(i,:)=X(i-1,:) + DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    F = Function('F', {X0}, {X}, {'x0'}, {'x'});
 
 foo = F('x0', [1;0]);