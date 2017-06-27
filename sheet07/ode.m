function dx = ode(t,x,u)

    dx = zeros(2,1);
    dx = [  x(2); ...
            sin(x(1)) + u(1)];
        
end

