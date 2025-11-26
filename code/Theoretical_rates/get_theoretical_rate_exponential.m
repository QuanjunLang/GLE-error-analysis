function rate = get_theoretical_rate_exponential(beta, gamma, C)


% rate = 0.5*((beta - gamma) + sqrt((beta + gamma)^2 -4*C));

rate = -0.5*((beta + gamma) - sqrt((beta - gamma).^2 + 4*C));
end