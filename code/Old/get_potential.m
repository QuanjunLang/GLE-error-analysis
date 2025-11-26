%% 
function V = get_potential(type, m)
    switch type
        case 'quadratic'
            % m = 20;
            V = @(x) m * (x - 5);
        case 'doublewell'
            V = @(x) 4 * x .* (x.^2 - 4);
        case 'cubic'
            V = @(x) 1 * (x - 5)^3;
    end
end


