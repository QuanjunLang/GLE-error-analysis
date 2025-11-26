close all
clear all
clc
addpaths


rng(1);
% dX_t = V_t * dt
% dV_t = -gamma * V_t * dt ...
%        - u * grad_U(X_t) * dt ...
%        - (integral from 0 to t of K(t, s) * V_s ds) * dt ...
%        + sigma * dB_t

%% System settings
sysInfo.d       = 3;
sysInfo.gamma   = 3;         % Damping coefficient
sysInfo.dt      = 0.02;          % Time step
sysInfo.T       = 50;             % Total simulation time
sysInfo.sigma   = 0;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);
%% Generate two kernels to compare
decay = 0.5;
amp1 = 1;
amp2 = 1.1;
intercept1 = 1;
intercept2 = 1;


kernel_type = 'exp';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions

plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log





%% Potential settings
% potential settings
potentialInfo.type          = 'quadratic_new';
potentialInfo.x_min_temp    = zeros(sysInfo.d, 1); % This is a fake x_min G where U = xRx + G. Will need to fit the true min
potentialInfo.u             = 1;    % u is just a parameter to adjust the Lyapunov function
potentialInfo.kappa         = 20;    % kappa is the minimal eigenvalue of the Hessian of the Potential

potentialInfo   = potential_info(sysInfo, potentialInfo);



%% Lyapunov settings
LG = potentialInfo.LG;
u = potentialInfo.u;
kappa = potentialInfo.kappa;

gamma = sysInfo.gamma;

p1 = LG * u * gamma^(-2);
fprintf('Requirement equation (4.7), we need L_G u gamma^(-2) < 3/4, where the actual parameter is %.4f\n', p1);



% lambda = 0.001;


p2 = kappa * u * gamma^(-2)/2;
lambda = min(p2, 1/8);
fprintf('Requirement equation (4.4), lambda < min(1/8, kappa u gamma^(-2)/2), where the actual parameter is %.6f and lambda is %.4f \n', p2, lambda);



s = 4;
mu = -2*decay+1;
% decay = 3;

RHS = 2*(2 * F_s_a(s, 2 * decay + mu ) / (s-1))^(1/2)
LHS = mu + 2 * gamma * lambda


%%


LyapunovInfo = Lyapunov_info(sysInfo, potentialInfo, lambda);



%%
get_theoretical_rate_exponential(-decay, 2*gamma*lambda, amp1)



%% Generate trajectory
noise = get_noise(sysInfo);

traj_input{1}.x0 = randn(sysInfo.d, 1);
traj_input{1}.v0 = randn(sysInfo.d, 1);
traj_input{1}.noise = noise;

trajInfo{1} = generate_data_second_order_highdim(sysInfo, kernelInfo.K1, potentialInfo, traj_input{1});

traj_input{2}.x0 = randn(sysInfo.d, 1);
traj_input{2}.v0 = randn(sysInfo.d, 1);
% 
% traj_input{2}.x0 = traj_input{1}.x0;
% traj_input{2}.v0 = traj_input{1}.v0;

traj_input{2}.noise = noise;
trajInfo{2} = generate_data_second_order_highdim(sysInfo, kernelInfo.K2, potentialInfo, traj_input{2});


%%
plot_trajectory_highdim(trajInfo{1}, sysInfo)
plot_trajectory_highdim(trajInfo{2}, sysInfo)


%%
traj_diff_Info.x = trajInfo{1}.x - trajInfo{2}.x;
traj_diff_Info.v = trajInfo{1}.v - trajInfo{2}.v;

plot_trajectory_decay_fit_new(traj_diff_Info, sysInfo, potentialInfo, 'trajectory difference', 1);


%%

r_square_1 = compute_Lyapunov_distance(trajInfo{1}.x, trajInfo{1}.v, potentialInfo.U_argmin, zeros(sysInfo.d, 1), LyapunovInfo);
r_square_2 = compute_Lyapunov_distance(trajInfo{2}.x, trajInfo{2}.v, potentialInfo.U_argmin, zeros(sysInfo.d, 1), LyapunovInfo);
r_square_diff = compute_Lyapunov_distance(trajInfo{1}.x, trajInfo{1}.v, trajInfo{2}.x, trajInfo{2}.v, LyapunovInfo);


fitInfo.t0 = 10;
fitInfo.v0 = sysInfo.sigma;
fit_new(r_square_diff, sysInfo, fitInfo);


%%
s = 4;
mu = -0.5;
% decay = 3;

RHS = 2*(2 * F_s_a(s, 2 * decay + mu ) / (s-1))^(1/2)
LHS = mu + 2 * gamma * lambda




% figure;
% sgrid = linspace(0, 10, 1000);
% plot(sgrid,  2*(2 * F_s_a(sgrid, 2 * decay + mu ) ./ (sgrid-1)).^(1/2))