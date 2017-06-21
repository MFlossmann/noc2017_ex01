function [ x_dot ] = J_FAD( u, m ,x0, nx , h )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x_dot = forw_AD(u, x0, nx, h);

tmp = size(x_dot,1) - m;
x_dot = x_dot(tmp:end , :);

end

