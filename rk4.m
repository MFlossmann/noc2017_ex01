function [ x_next ] = rk4( h, x0, u, ode, param )
%rk4 Computes one Runge-Kutta-method (RK4) step
%   Detailed explanation goes here
k1 = ode(h, x0, u, param);
k2 = ode(h, x0 + k1 * h/2, u, param);
k3 = ode(h, x0 + k2 * h/2, u, param);
k4 = ode(h, x0 + k3 * h, u, param);

x_next = x0 + h * (k1 + 2*k2 + 2*k3 + k4) / 6;

end

