omega = 10;
beta = 8/3;
phi = 28;

x_0 = [1; 0; 0];                            % initial condition
h = 0.01;                                   % step size
t = [0:h:100; 0:h:100; 0:h:100];            % 
x = zeros(3, length(t));
x(:,1) = x_0;

%Declaring the function
F = @(x)[omega*(x(2)-x(1)); x(1)*(phi-x(3))-x(2); x(1)*x(2)-beta*x(3)]; 

%RK 4
for i=1 : (length(t)-1)
    k_1 = F(x(:,i));
    k_2 = F(x(:,i)+0.5*h*k_1);
    k_3 = F(x(:,i)+0.5*h*k_2);
    k_4 = F(x(:,i)+h*k_3);
    
    x(:,i+1) = x(:,i)+(h/6)*(k_1+2*k_2+2*k_3+k_4);
end

