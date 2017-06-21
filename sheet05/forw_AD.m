function [ x_dot ] = forw_AD( u ,x0, nx , h)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

N = size(u,1) + 1;

% initialize w and w_dot
w = zeros(7,1);
w_dot = zeros(7 ,N-1);

x = [x0];
x_dot = zeros(N-1,N-1);

for i=1:N-1
    % first define the subfunctions
    w(1) = x(end);
    w(2) = u(i);
    w(3) = 1 - w(1);
    w(4) = w(3) * w(1);
    w(5) = w(4) + w(2);
    w(6) = h*w(5);
    w(7) = w(1) + w(6);
    
    % append w(7) as most recent x
    x = [x, w(7)];
    
    w_dot(1,:) = x_dot(i,:);
    w_dot(2,:) = zeros(1,N-1); w_dot(2,i) = 1;
    w_dot(3,:) = -w_dot(1,:);
    w_dot(4,:) = w(1)*w_dot(3,:) + w_dot(1,:)*w(3);
    w_dot(5,:) = w_dot(4,:)+ w_dot(2,:);
    w_dot(6,:) = h*w_dot(5,:);
    w_dot(7,:) = w_dot(1,:) + w_dot(6,:);
    
    x_dot(i+1,:) = w_dot(7,:);
end

end