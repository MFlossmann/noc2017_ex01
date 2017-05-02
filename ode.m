function [ xdot ] = ode( t, x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xdot(1) = x(2);
xdot(2) = -0.2*x(2) - x(1);

end

