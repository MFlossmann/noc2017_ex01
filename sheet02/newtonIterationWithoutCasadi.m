clear; clc; close;

iteration_length = 1000;
initial_value = [1; 1.1];
phi = 500;
approx_hess = [1/phi 0;0 1/phi];

w = zeros(2,iteration_length);
v = zeros(2,iteration_length);
w(:,1) = initial_value;
v(:,1) = initial_value;

for i=1:1:iteration_length
    [grad_w,hess_w] = gradientHessianFunction(w(1,i),w(2,i));
    w(:,i+1)=w(:,i)-hess_w^(-1)*grad_w;
end

for i=1:1:iteration_length
    grad_v = gradientHessianFunction(v(1,i),v(2,i));
    v(:,i+1)= v(:,i)-approx_hess*grad_v;
end

figure;
hold on;
plot(v(:,1),v(:,2), '--');
plot(w(:,1),w(:,2), ':');
hold off;
