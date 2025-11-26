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
sysInfo.d       = 2;
sysInfo.gamma   = 10;         % Damping coefficient
sysInfo.dt      = 0.05;          % Time step
sysInfo.T       = 300;             % Total simulation time
sysInfo.sigma   = 0;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);
%% Generate two kernels to compare
decay = 0.1;
amp1 = 1;
amp2 = 1.1;
intercept1 = 0.01;
intercept2 = 0.01;


kernel_type = 'exp';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions

plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log





%% Potential settings
potential_type = 'quadratic_new';
u = 3;


potentialInfo   = potential_info(potential_type, sysInfo, 'x_min', [0, 0]', 'u', u);
plot_potential_2d(potentialInfo, [-5, 5], [-5, 5], 200);


%% Lyapunov settings
LG = potentialInfo.LG;
u = potentialInfo.u;
kappa = potentialInfo.kappa;

gamma = sysInfo.gamma;

p1 = LG * u * gamma^(-2);
fprintf('Requirement equation (4.7), we need L_G u gamma^(-2) < 3/4, where the actual parameter is %.4f\n', p1);



lambda = 0.001;


p2 = min(kappa * u * gamma^(-2)/2, 1/8);
fprintf('Requirement equation (4.4), lambda < min(1/8, kappa u gamma^(-2)/2), where the actual parameter is %.6f and lambda is %.4f \n', p2, lambda);
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

plot_trajectory_decay_fit_new(traj_diff_Info, sysInfo, potentialInfo, 'trajectory difference', 1)


%%

r_square_1 = compute_Lyapunov_distance(trajInfo{1}.x, trajInfo{1}.v, potentialInfo.U_argmin, [0, 0]', LyapunovInfo);
r_square_2 = compute_Lyapunov_distance(trajInfo{2}.x, trajInfo{2}.v, potentialInfo.U_argmin, [0, 0]', LyapunovInfo);
r_square_diff = compute_Lyapunov_distance(trajInfo{1}.x, trajInfo{1}.v, trajInfo{2}.x, trajInfo{2}.v, LyapunovInfo);


fitInfo.t0 = 10;
fitInfo.v0 = 1e-100;
fit_new(r_square_diff, sysInfo, fitInfo)
%%




% 
% 
% 
% % Time vector (assuming consistent time steps)
% T = length(r_square_1);
% t = sysInfo.t_grid;
% 
% % Plot
% figure;
% 
% % 1. Regular scale
% subplot(1,3,1); hold on;
% plot(t, r_square_1, 'b-', 'DisplayName', '$r((X_t, V_t), (x^*, v^*))^2$');
% plot(t, r_square_2, 'r--', 'DisplayName', '$r((\tilde X_t, \tilde V_t), (x^*, v^*))^2$');
% plot(t, r_square_diff, 'k-.', 'DisplayName', '$r((X_t, V_t), (\tilde X_t, \tilde V_t))^2$');
% xlabel('Time step'); ylabel('Lyapunov value');
% title('Regular scale');
% legend('Interpreter','latex')
% set(gca, 'FontSize', 18)
% grid on;
% 
% 
% 
% % 2. Semilogy
% subplot(1,3,2);hold on;
% semilogy(t, r_square_1, 'b-', 'DisplayName', '$r((X_t, V_t), (x^*, v^*))^2$');
% semilogy(t, r_square_2, 'r--', 'DisplayName', '$r((\tilde X_t, \tilde V_t), (x^*, v^*))^2$');
% semilogy(t, r_square_diff, 'k-.', 'DisplayName', '$r((X_t, V_t), (\tilde X_t, \tilde V_t))^2$');
% xlabel('Time step'); ylabel('Lyapunov value');
% title('Semilogy');
% legend('Interpreter','latex')
% set(gca, 'FontSize', 18, 'YScale', 'log')
% grid on;
% 
% 
% 
% % 3. Log-log
% subplot(1,3,3);hold on;
% loglog(t, r_square_1, 'b-', 'DisplayName', '$r((X_t, V_t), (x^*, v^*))^2$');
% loglog(t, r_square_2, 'r--', 'DisplayName', '$r((\tilde X_t, \tilde V_t), (x^*, v^*))^2$');
% loglog(t, r_square_diff, 'k-.', 'DisplayName', '$r((X_t, V_t), (\tilde X_t, \tilde V_t))^2$');
% xlabel('Time step'); ylabel('Lyapunov value');
% title('Log-log'); 
% legend('Interpreter','latex')
% set(gca, 'FontSize', 18, 'YScale', 'log', 'XScale', 'log')
% grid on;
% 
% 
% 
% %%
% % Extract trajectories
% x1 = trajInfo{1}.x;  % [d x N]
% v1 = trajInfo{1}.v;
% x2 = trajInfo{2}.x;
% v2 = trajInfo{2}.v;
% 
% % Time vector
% T = size(x1, 2);
% t = sysInfo.t_grid;
% 
% % Compute component-wise squared distances
% x_diff_sq = sum((x1 - x2).^2, 1);     % Position difference squared
% v_diff_sq = sum((v1 - v2).^2, 1);     % Velocity difference squared
% traj_diff_sq = x_diff_sq + v_diff_sq; % Lyapunov-like total difference
% 
% % Plotting
% figure;
% 
% % ----- (1) Regular Scale -----
% subplot(1, 3, 1);
% plot(t, x_diff_sq, 'b-', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta x\|^2$');
% hold on;
% plot(t, v_diff_sq, 'r--', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta v\|^2$');
% plot(t, traj_diff_sq, 'k-.', 'LineWidth', 2, 'DisplayName', '$\|\Delta x\|^2 + \|\Delta v\|^2$');
% xlabel('Time $t$', 'Interpreter', 'latex');
% ylabel('Squared Distance', 'Interpreter', 'latex');
% title('Regular Scale', 'Interpreter', 'latex');
% legend('Interpreter', 'latex', 'Location', 'best'); grid on;
% 
% % ----- (2) Semilog-y Scale -----
% subplot(1, 3, 2);
% semilogy(t, x_diff_sq, 'b-', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta x\|^2$');
% hold on;
% semilogy(t, v_diff_sq, 'r--', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta v\|^2$');
% semilogy(t, traj_diff_sq, 'k-.', 'LineWidth', 2, 'DisplayName', '$\|\Delta x\|^2 + \|\Delta v\|^2$');
% xlabel('Time $t$', 'Interpreter', 'latex');
% ylabel('Squared Distance (log scale)', 'Interpreter', 'latex');
% title('Semilog-y Scale', 'Interpreter', 'latex');
% legend('Interpreter', 'latex', 'Location', 'best'); grid on;
% 
% % ----- (3) Log-Log Scale -----
% subplot(1, 3, 3);
% loglog(t, x_diff_sq, 'b-', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta x\|^2$');
% hold on;
% loglog(t, v_diff_sq, 'r--', 'LineWidth', 1.5, 'DisplayName', '$\|\Delta v\|^2$');
% loglog(t, traj_diff_sq, 'k-.', 'LineWidth', 2, 'DisplayName', '$\|\Delta x\|^2 + \|\Delta v\|^2$');
% xlabel('Time $t$', 'Interpreter', 'latex');
% ylabel('Squared Distance (log-log)', 'Interpreter', 'latex');
% title('Log-Log Scale', 'Interpreter', 'latex');
% legend('Interpreter', 'latex', 'Location', 'best'); grid on;
% 
% % Optional: set global figure properties
% sgtitle('Trajectory Difference Analysis (Lyapunov vs Euclidean)', 'Interpreter', 'latex');
