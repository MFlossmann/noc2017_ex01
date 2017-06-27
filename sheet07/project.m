function index = project( in, values )
    % ASSUMPTION: the argument "values" contains a vector of equally spaced points
    min = values(1);
    max = values(end);
    N = length(values);
    
    index = round(((in-min)/(max-min))*(N-1)) + 1;
end

