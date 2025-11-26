close all
clear all
clc
addpaths
rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 2;
sysInfo.gamma   = 1;             % Damping coefficient
sysInfo.dt      = 0.02;          % Time step
sysInfo.T       = 200;           % Total simulation time

sysInfo.sigma   = 1e-6;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);

%% Generate two kernels to compare
kernel_type = 'osc';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
decay = 2;
amp1 = 1;
amp2 = 1.1;
intercept1 = 0.1;
intercept2 = 0.1;

kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions
plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log




rate = get_theoretical_rate_exponential(decay, sysInfo.gamma, amp1);
fprintf('\nFor exponential kernel, the theoretical rate is : %.2f\n', rate)

%% Generate trajectory
noise = get_noise(sysInfo);

traj_input{1}.v0 = randn(sysInfo.d, 1);
traj_input{1}.noise = noise;



trajInfo{1} = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});



traj_input{2}.v0 = randn(sysInfo.d, 1);
traj_input{2}.noise = noise;
trajInfo{2} = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});


%% 
plot_velocity_highdim(trajInfo{1}.v, sysInfo, 'traj 1', true);
plot_velocity_highdim(trajInfo{2}.v, sysInfo, 'traj 2', true);



fitInfo = plot_velocity_highdim(abs(trajInfo{1}.v - trajInfo{2}.v), sysInfo, 'traj difference', true);
%%


% 
% 
% 
% 
% % Parameters
% C = 1;                 % example value
% intercept = 0.1;         % example value
% beta_vals = linspace(0.01, 3, 500);  % Avoid beta = 1 to prevent singularity
% 
% % Compute thresholds
% threshold_vals = arrayfun(@(b) get_theoretical_threshold_powerlaw(b, C, intercept), beta_vals);
% 
% % Plot
% figure;
% plot(beta_vals, threshold_vals, 'LineWidth', 2);
% xlabel('\beta');
% ylabel('Threshold');
% title('Theoretical Threshold vs. \beta');
% grid on;
% xline(1, '--r', 'Singularity at \beta = 1');
% ylim([-10, 50])