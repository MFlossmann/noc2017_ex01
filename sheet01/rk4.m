function [ x_sim ] = rk4(T, h, ode, x0)
%rk4 Computes one Runge-Kutta-method (RK4) step
%   Detailed explanation goes here

N = T/h;

x_sim = zeros(size(x0,1),N);

x_sim(:,1) = x0;

for i = 2:N
    k1 = ode(h*i, x_sim(:,i-1));
    k2 = ode(h*i, x_sim(:,i-1) + k1 * h/2);
    k3 = ode(h*i, x_sim(:,i-1) + k2 * h/2);
    k4 = ode(h*i, x_sim(:,i-1) + k3 * h);

    x_sim(:,i) = x_sim(:,i-1) + h * (k1 + 2*k2 + 2*k3 + k4) / 6;
end

x_sim = x_sim';

end

