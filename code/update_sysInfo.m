function sysInfo = update_sysInfo(sysInfo)


% sigma = sysInfo.sigma;
% dt = sysInfo.dt;
d = sysInfo.d;


sysInfo.N       = sysInfo.T / sysInfo.dt;   % Number of time steps


sysInfo.t_grid       = (0:sysInfo.N-1) * sysInfo.dt;  % Time vector
% sysInfo.noise   = sqrt(sigma^2 / dt) * randn(1, sysInfo.N); % Shared Brownian noise

sysInfo.x0 = randn(d, 1);
sysInfo.v0 = randn(d, 1);
sysInfo.tn = length(sysInfo.t_grid);
end

