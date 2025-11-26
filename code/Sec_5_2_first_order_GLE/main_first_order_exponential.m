close all
clear all
clc
addpaths
rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 1;
sysInfo.gamma   = 3;         % Damping coefficient
sysInfo.dt      = 0.02;          % Time step
sysInfo.T       = 30;             % Total simulation time

sysInfo.sigma   = 1e-4;

sysInfo         = update_sysInfo(sysInfo);

%% Generate two kernels to compare
kernel_type = 'osc';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
decay = 3;
amp1 = 2;
amp2 = 2.1;
intercept1 = 0.1;
intercept2 = 0.1;

kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions
plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log
y



rate = get_theoretical_rate_exponential(decay, sysInfo.gamma, amp1);
fprintf('\nFor exponential kernel, the theoretical rate is : %.2f\n', rate)

%% Generate Ensemble trajectories
M = 10;
traj_1 = zeros(M, sysInfo.tn);
traj_2 = zeros(M, sysInfo.tn);

v0_sigma = 0.5;
v0_mu = 1;

for i = 1:M

    noise = get_noise(sysInfo);

    traj_input{1}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{1}.noise = noise;

    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});
    traj_1(i, :) = temp.v;

    traj_input{2}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{2}.noise = noise;
    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});
    traj_2(i, :) = temp.v;
end

%%

traj_1_mean = mean(traj_1.^2);
traj_2_mean = mean(traj_2.^2);
traj_diff = mean((traj_1 - traj_2).^2);
c2 = min(traj_diff);
%%
fitOptions.fit_t0 = 2;
fitOptions.fit_v0 = 1e-8;
fitOptions.plotON = true;
fitOptions.type = 'semilog';    % semilog (Exp) or loglog (Powerlaw) or all


fitOptions_1        = fitOptions;
fitOptions_1.ttl    = 'traj 1';

fitInfo_1 = fit_traj_rate(traj_1_mean, sysInfo, fitOptions_1);

fitOptions_2        = fitOptions;
fitOptions_2.ttl    = 'traj 2';
fitInfo_2 = fit_traj_rate(traj_2_mean, sysInfo, fitOptions_2);

fitOptions_diff        = fitOptions;
fitOptions_diff.ttl    = 'traj difference';

fitInfo_diff = fit_traj_rate(traj_diff, sysInfo, fitOptions_diff);


%%

% get_theoretical_rate(decay, sysInfo.gamma, amp1)





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