function [ x_next ] = euler_step( h, x0, u, ode, param )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% M = Magic number... = 1 for this exercise
M = 1;

h = h / M;

x_next = x0;
for i=1:M
    x_next = x_next + h * ode(h, x_next);
end

end
