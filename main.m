%% 2. Explicit integrators

%% 2.a
x0 = [1;0]; % this is a column vector
T_s = 0.5;
k = 0:1:19; % to create an array from 0 to 19
t = k*T_s; % This means that the scalar T_s is multiplied with all entries of k

%[t,x] = ode45(@(t,x) [x(2); -0.2*x(2) - x(1)] , t, x0);
[t,x] = ode45(ode, t, x0);

figure;
hold on;
plot(x(:,1),x(:,2),'-o');
hold off;

%% 2.b

ode = diff(x,t) == [x(2); -0.2*x(2) - x(1)];