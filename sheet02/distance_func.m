function [ distance ] = distance_func( v_0 )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

P_t = ballistic_dynamics_RK4(v_0);

distance = norm(P_t(1:2) - P_t(3:4));

end

