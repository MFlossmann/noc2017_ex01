%%Task 1.c
clear; clc; close;

iteration_length = 1000;
initial_value = [1; 1.1];
rho = 500;
approx_hess = eye(2).*rho;
approx_inv = inv(approx_hess);

% w: analytic hessian
% v: approximated hessian
w = zeros(2,iteration_length);
v = zeros(2,iteration_length);
w(:,1) = initial_value;
v(:,1) = initial_value;

for i=1:1:iteration_length
    w(:,i+1)= w(:,i) - inv(rosenbrock_hessian(w(1,i),w(2,i)))...
                       * rosenbrock_grad(w(1,i),w(2,i));
end

for i=1:1:iteration_length
    grad_v = gradientHessianFunction(v(1,i),v(2,i));
    v(:,i+1)= v(:,i) - approx_inv*rosenbrock_grad(v(1,i),v(2,i));
end

figure;
hold on;
plot(v(1,:),v(2,:), 'rx-');
plot(w(2,:),w(2,:), 'b+-');
legend('M=\rho I','M=\nabla^2 f(x,y)')legend('M=\rho I','M=\nabla^2 f(x,y)')

hold off;

%% Task 1.d
clear all;
import casadi.*

x = MX.sym('x');
y = MX.sym('y');

f = Function('f', {x,y}, {(1-x)^2 + 100*(y - x^2)^2});

F = (1-x)^2 + 100*(y - x^2)^2;

j = jacobian(F,[x,y]);
h = hessian(F,[x,y]);