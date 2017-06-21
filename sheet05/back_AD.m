function [ x_dot ] = back_AD( u ,x0, nx , h )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = size(u,1);

w = zeros(7,1);
w_bar = zeros(7, N);

x_dot = zeros(N,1);

% init seed vector
w_bar(end) = 1;

x = x0;

for i=1:N
    % define functions
    w(1) = x;
    w(2) = u(i);
    w(3) = 1 - w(1);
    w(4) = w(3) * w(1);
    w(5) = w(4) + w(2);
    w(6) = h*w(5);
    w(7) = w(1) + w(6);
    
    x = w(7);
    
    % backwards sweep
    % 7
    w_bar(1) = w_bar(1) + w_bar(7);
    w_bar(6) = w_bar(6) + w_bar(7);
    
    % 6
    w_bar(5) = w_bar(5) + h*w_bar(6);
    
    % 5
    w_bar(4) = w_bar(4) + 1*w_bar(5);
    w_bar(2) = w_bar(2) + w_bar(5);
    
    % 4
    w_bar(3) = w_bar(3) + w(1) * w_bar(4);
    w_bar(1) = w_bar(1) + w(3) * w_bar(4);
    
    % 3
    w_bar(1) = w_bar(1) - w_bar(3);
    
end

end

