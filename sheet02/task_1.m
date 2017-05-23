%%Task 1.c
clear all; clc; close;

iteration_length = 1000;
initial_value = [1; 1.1];
rho = 500;
approx_hess = eye(2).*rho;
approx_inv = inv(approx_hess);

% w: analytic hessian
% v: approximated hessian
w = zeros(2,iteration_length);
v = zeros(2,iteration_length);
w(:,1) = initial_value(:);
v(:,1) = initial_value(:);

for i=1:1:iteration_length
    w(:,i+1)= w(:,i) - inv(rosenbrock_hessian(w(1,i),w(2,i)))...
                       * rosenbrock_grad(w(1,i),w(2,i));
end

for i=1:1:iteration_length
    v(:,i+1)= v(:,i) - approx_inv*rosenbrock_grad(v(1,i),v(2,i));
end

fig1 = figure;
hold on;
plot(v(1,:),v(2,:), 'rx-');
plot(w(1,:),w(2,:), 'b+-');
legend('M=\rho I','M=\nabla^2 f(x,y)')
hold off;

%% Task 1.d
import casadi.*

x = MX.sym('x',2);

%f = Function('f', {x,y}, {(1-x)^2 + 100*(y - x^2)^2}); 

% Rosenbrock function
F = (1-x(1))^2 + 100*(x(2) - x(1)^2)^2;

grad = Function('grad', {x}, {jacobian(F,x)});
H = Function('H', {x}, {hessian(F,x)});

u0 = DM([1; 1.1]);
u = DM(2,iteration_length);

u(:,1) = u0;

for i=2:iteration_length
    u(:,i) = u(:,i-1) - inv(H(u(:,i-1))) *...
                        grad(u(:,i-1))';
end

trace = full(u);

fig1;
hold on;
plot(trace(1,:),trace(2,:), 'g.-');
hold off;