function noise = get_noise(sysInfo)

sigma = sysInfo.sigma;
dt = sysInfo.dt;


noise = sqrt(sigma^2 / dt) * randn(1, sysInfo.N);

end