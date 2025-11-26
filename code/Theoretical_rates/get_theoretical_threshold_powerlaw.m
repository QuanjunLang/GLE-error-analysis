function threshold = get_theoretical_threshold_powerlaw(beta, C, intercept)

threshold = C/(beta-1) * intercept ^(1- beta);


end