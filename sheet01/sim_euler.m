function [ x_sim ] = sim_euler( t, h, xdot, x0 )
%sim_euler Simulates the path of a robot, given the initial position and
%control variable u

N = t/h;
m = size(x0,1);

x_sim = zeros(m,N);
x_sim(:,1) = x0;

for i = 2:N
    x_sim(:,i) = x_sim(:,i-1) + xdot(i*h, x_sim(:,i-1));
end

x_sim = x_sim';

end
