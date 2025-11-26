function sysInfo = system_info()
% System parameters

sysInfo.mass = 1;
sysInfo.gamma = 10;         % Damping coefficient
sysInfo.dt = 0.01;          % Time step
sysInfo.T = 50;             % Total simulation time
sysInfo.N = sysInfo.T / sysInfo.dt;   % Number of time steps
sysInfo.t = (0:sysInfo.N-1) * sysInfo.dt;  % Time vector

end