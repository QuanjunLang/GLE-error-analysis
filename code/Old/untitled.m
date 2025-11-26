close all; 
clear all; clc;

%% Parameters
mass = 1;
gamma = 10;
dt = 0.04;
T = 12;
N = round(T / dt);
t = (0:N-1) * dt;

%% Choose memory kernels
kernel_type = 'power';  % options: 'exp', 'power', 'osc'
[f1, f2, K1, K2] = switch_kernel(kernel_type);

%% Choose external potential
V_type = 'quadratic';  % options: 'quadratic', 'doublewell', 'linear', 'sinusoidal'
m = 20;
V = get_potential(V_type, m);

%% Shared noise
sigma = 0.001;
eta = sqrt(2 * gamma * sigma^2 / dt) * randn(1, N);

%% Initial conditions
x1_0 = 6; v1_0 = 2;
x2_0 = 4; v2_0 = -2;

x1 = zeros(1, N); v1 = zeros(1, N);
x2 = zeros(1, N); v2 = zeros(1, N);
x1(1) = x1_0; v1(1) = v1_0;
x2(1) = x2_0; v2(1) = v2_0;

% %% Memory integral storage
% mem1_vals = zeros(1, N);
% mem2_vals = zeros(1, N);

%% Time stepping with RK4
for n = 2:N
    tn = t(n);

    % % Compute memory integrals
    % if n > 2
    %     t_hist = t(1:n-1);
    %     v1_hist = v1(n-1:-1:1);
    %     v2_hist = v2(n-1:-1:1);
    % 
    %     K1_vals = arrayfun(@(s) K1(tn, s), t_hist);
    %     K2_vals = arrayfun(@(s) K2(tn, s), t_hist);
    % 
    %     mem1_vals(n) = trapz(t_hist, K1_vals .* v1_hist);
    %     mem2_vals(n) = trapz(t_hist, K2_vals .* v2_hist);
    % else
    %     mem1_vals(n) = 0;
    %     mem2_vals(n) = 0;
    % end

    if n == 2
        % Use Euler step to start
        mem1 = 0;
        mem2 = 0;
        v1(n) = v1(n-1) + dt/mass * (-gamma * v1(n-1) - V(x1(n-1)) - mem1 + eta(n));
        v2(n) = v2(n-1) + dt/mass * (-gamma * v2(n-1) - V(x2(n-1)) - mem2 + eta(n));
        x1(n) = x1(n-1) + v1(n) * dt;
        x2(n) = x2(n-1) + v2(n) * dt;
        continue;
    end

    %% RK4 for trajectory 1
    kx1 = v1(n-1);
    kv1 = RHS(x1(n-1), v1(n-1), tn, v1, t, K1, V, dt, n, gamma);

    kx2 = v1(n-1) + 0.5 * dt * kv1;
    kv2 = RHS(x1(n-1) + 0.5 * dt * kx1, v1(n-1) + 0.5 * dt * kv1, tn, v1, t, K1, V, dt, n, gamma);

    kx3 = v1(n-1) + 0.5 * dt * kv2;
    kv3 = RHS(x1(n-1) + 0.5 * dt * kx2, v1(n-1) + 0.5 * dt * kv2, tn, v1, t, K1, V, dt, n, gamma);

    kx4 = v1(n-1) + dt * kv3;
    kv4 = RHS(x1(n-1) + dt * kx3, v1(n-1) + dt * kv3, tn, v1, t, K1, V, dt, n, gamma);

    x1(n) = x1(n-1) + (dt / 6) * (kx1 + 2*kx2 + 2*kx3 + kx4);
    v1(n) = v1(n-1) + (dt / 6) * (kv1 + 2*kv2 + 2*kv3 + kv4) + eta(n);

    %% RK4 for trajectory 2
    kx1 = v2(n-1);
    kv1 = RHS(x2(n-1), v2(n-1), tn, v2, t, K2, V, dt, n, gamma);

    kx2 = v2(n-1) + 0.5 * dt * kv1;
    kv2 = RHS(x2(n-1) + 0.5 * dt * kx1, v2(n-1) + 0.5 * dt * kv1, tn, v2, t, K2, V, dt, n, gamma);

    kx3 = v2(n-1) + 0.5 * dt * kv2;
    kv3 = RHS(x2(n-1) + 0.5 * dt * kx2, v2(n-1) + 0.5 * dt * kv2, tn, v2, t, K2, V, dt, n, gamma);

    kx4 = v2(n-1) + dt * kv3;
    kv4 = RHS(x2(n-1) + dt * kx3, v2(n-1) + dt * kv3, tn, v2, t, K2, V, dt, n, gamma);

    x2(n) = x2(n-1) + (dt / 6) * (kx1 + 2*kx2 + 2*kx3 + kx4);
    v2(n) = v2(n-1) + (dt / 6) * (kv1 + 2*kv2 + 2*kv3 + kv4) + eta(n);
end

%% Compute |x1 - x2|
x_diff = abs(x1 - x2);

%% Use tiledlayout for 4 subplots
figure('Position', [100, 100, 1000, 800]);
tl = tiledlayout(4,1, 'TileSpacing', 'compact');

% 1. Position plot
nexttile;
plot(t, x1, 'b-', 'LineWidth', 2, 'DisplayName', 'x_1(t)');
hold on;
plot(t, x2, 'r--', 'LineWidth', 2, 'DisplayName', 'x_2(t)');
ylabel('Position');
title('Position Trajectories');
legend('Location', 'best');
grid on;

% 2. Velocity plot
nexttile;
plot(t, v1, 'b-.', 'LineWidth', 2, 'DisplayName', 'v_1(t)');
hold on;
plot(t, v2, 'r:', 'LineWidth', 2, 'DisplayName', 'v_2(t)');
ylabel('Velocity');
title('Velocity Trajectories');
legend('Location', 'best');
grid on;

% 3. Difference |x1 - x2|
nexttile;
plot(t, abs(x1 - x2), 'k', 'LineWidth', 2);
ylabel('|x_1 - x_2|');
title('Trajectory Difference');
% set(gca, 'YScale', 'log');
grid on;

% 4. Memory kernels
nexttile;
tau = linspace(0, T, 1000);
plot(tau, f1(tau), 'b', 'LineWidth', 2, 'DisplayName', 'f_1(\tau)');
hold on;
plot(tau, f2(tau), 'r--', 'LineWidth', 2, 'DisplayName', 'f_2(\tau)');
xlabel('\tau'); ylabel('Kernel Value');
title('Memory Kernels');
legend('Location', 'best');
grid on;


%%
% Absolute-value log-scale plots for position and velocity
figure('Position', [100, 100, 1000, 800]);
tiledlayout(4,1, 'TileSpacing', 'compact');

% 1. |x1| and |x2|
nexttile;
plot(t, abs(x1), 'b', 'LineWidth', 2, 'DisplayName', '|x_1(t)|');
hold on;
plot(t, abs(x2), 'r--', 'LineWidth', 2, 'DisplayName', '|x_2(t)|');
yscale('log');
ylabel('|x(t)|');
title('Log-Scale Position Magnitudes');
legend('Location', 'best');
grid on;

% 2. |v1| and |v2|
nexttile;
plot(t, abs(v1), 'b-.', 'LineWidth', 2, 'DisplayName', '|v_1(t)|');
hold on;
plot(t, abs(v2), 'r:', 'LineWidth', 2, 'DisplayName', '|v_2(t)|');
yscale('log');
ylabel('|v(t)|');
title('Log-Scale Velocity Magnitudes');
legend('Location', 'best');
grid on;

% 3. |x1 - x2|
nexttile;
plot(t, abs(x1 - x2), 'k', 'LineWidth', 2);
yscale('log');
ylabel('|x_1 - x_2|');
title('Log-Scale Position Difference');
grid on;

% 4. |v1 - v2|
nexttile;
plot(t, abs(v1 - v2), 'k--', 'LineWidth', 2);
yscale('log');
ylabel('|v_1 - v_2|');
xlabel('Time');
title('Log-Scale Velocity Difference');
grid on;

%%
% Lyapunov-like terms
lambda = min(0.125, m / (gamma^2));

r1 = (m / gamma^2) * (x1 - x2).^2;
r2 = 0.5 * ((1 - 2 * lambda) * (x1 - x2) + (v1 - v2) / gamma).^2;
r3 = 0.5 / (gamma^2) * (v1 - v2).^2;
r = r1 + r2 + r3;

% Plot with formula-based legend
figure('Position', [100, 100, 900, 500]); hold on;

plot(t, log10(r1), 'b-',  'LineWidth', 2);
plot(t, log10(r2), 'r--', 'LineWidth', 2);
plot(t, log10(r3), 'k-.', 'LineWidth', 2);
plot(t, log10(r),  'm-',  'LineWidth', 2);

xlabel('Time'); ylabel('$\log_{10}(r)$', 'Interpreter', 'latex');
title('Lyapunov Function Components', 'Interpreter', 'latex');

lgd = legend({ ...
    '$r_1(t) = \frac{m}{\gamma^2}(x_1 - x_2)^2$', ...
    '$r_2(t) = \frac{1}{2}\left[(1 - 2\lambda)(x_1 - x_2) + \frac{1}{\gamma}(v_1 - v_2)\right]^2$', ...
    '$r_3(t) = \frac{1}{2\gamma}(v_1 - v_2)^2$', ...
    '$r(t) = r_1 + r_2 + r_3$' ...
}, 'Interpreter', 'latex', 'Location', 'best');

set(lgd, 'FontSize', 14); 
grid on;


%%




%% Memory integral storage
mem1_vals = zeros(1, N);
mem2_vals = zeros(1, N);

% Time stepping with RK4
for n = 2:N
    tn = t(n);

    % Compute memory integrals
    if n > 2
        t_hist = t(1:n-1);
        v1_hist = v1(1:1:n-1);
        v2_hist = v2(1:1:n-1);

        K1_vals = arrayfun(@(s) K1(tn, s), t_hist);
        K2_vals = arrayfun(@(s) K2(tn, s), t_hist);

        mem1_vals(n) = trapz(t_hist, K1_vals .* (v1_hist));
        mem2_vals(n) = trapz(t_hist, K2_vals .* (v2_hist));
    else
        mem1_vals(n) = 0;
        mem2_vals(n) = 0;
    end

end


figure;
plot(t, mem1_vals, 'b-', 'LineWidth', 2, 'DisplayName', '$\int_0^t K_1(t,s) v_1(s) ds$');
hold on;
plot(t, mem2_vals, 'r--', 'LineWidth', 2, 'DisplayName', '$\int_0^t K_2(t,s) v_2(s) ds$');
xlabel('Time', 'Interpreter', 'latex');
ylabel('Memory Integral', 'Interpreter', 'latex');
title('Memory Term Comparison', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;

%% 
z = x1 - x2;
w = v1 - v2;

% memory = gamma^(-1)*((1-2*lambda)*z + 2*gamma^(-1)*w).*(mem1_vals - mem2_vals);

% Compute time derivative of r(t)
dr_dt = gradient(r, dt);

% Reference decay expression
% rhs = -2 * lambda * gamma * r - gamma^(-1)*((1-2*lambda)*z + 2*gamma^(-1)*w).*(mem1_vals - mem2_vals);  % memory difference drives correction

lhs = gamma*(dr_dt + 2*lambda*gamma*r);


rhs = ((1-2*lambda)*z + 2*gamma^(-1)*w).*(mem1_vals - mem2_vals);


% rhs = ((1-2*lambda)*z + 2*gamma^(-1)*w).^2 + 0.5*mem1_vals.^2 + 0.5*mem2_vals.^2;


% Plot
figure('Position', [100, 100, 800, 500]); hold on;
plot(t, lhs, 'b-', 'LineWidth', 2);
plot(t, rhs,    'r--', 'LineWidth', 2);

xlabel('Time', 'Interpreter', 'latex');
ylabel('Value', 'Interpreter', 'latex');
title('Evolution of Lyapunov Function', 'Interpreter', 'latex');

legend({'lhs', 'rhs'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

grid on;


%% difference in memory

mem_diff_v1_vals = zeros(1, N);

% Time stepping with RK4
for n = 2:N
    tn = t(n);

    % Compute memory integrals
    if n > 2
        t_hist = t(1:n-1);
        v1_hist = v1(1:1:n-1);

        K_diff_vals = arrayfun(@(s) K1(tn, s) - K2(tn, s), t_hist);

        mem_diff_v1_vals(n) = trapz(t_hist, K_diff_vals .* (v1_hist));
    else
        mem_diff_v1_vals(n) = 0;
    end

end


mem2_v_diff_vals = zeros(1, N);

% Time stepping with RK4
for n = 2:N
    tn = t(n);

    % Compute memory integrals
    if n > 2
        t_hist = t(1:n-1);
        v_diff_hist = v1(1:1:n-1) - v2(1:1:n-1);

        K2_vals = arrayfun(@(s) K2(tn, s), t_hist);

        mem2_v_diff_vals(n) = trapz(t_hist, K2_vals .* (v_diff_hist));
    else
        mem2_v_diff_vals(n) = 0;
    end

end


lhs = gamma*(dr_dt + 2*lambda*gamma*r);

epsilon = 0.1;

rhs_new = epsilon*((1-2*lambda)*z + 2*gamma^(-1)*w).^2 + 1/epsilon*0.5*mem2_v_diff_vals.^2 + 1/epsilon*0.5*mem2_v_diff_vals.^2;



rhs_new_4r = 4*epsilon*r + 1/epsilon*0.5*mem2_v_diff_vals.^2 +  1/epsilon*0.5*mem2_v_diff_vals.^2;

% Plot
figure('Position', [100, 100, 800, 500]); hold on;
plot(t, lhs, 'b-', 'LineWidth', 2);
plot(t, rhs, 'k-', 'LineWidth', 2);
plot(t, rhs_new,    'r--', 'LineWidth', 2);
plot(t, rhs_new_4r,    'g--', 'LineWidth', 2);


xlabel('Time', 'Interpreter', 'latex');
ylabel('Value', 'Interpreter', 'latex');
title('Evolution of Lyapunov Function', 'Interpreter', 'latex');

legend({'lhs', 'rhs', 'rhs new', 'rhs new 4r'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

grid on;


figure;hold on;
plot(t, rhs);
plot(t, rhs_new);
set(gca, 'YScale', 'log');
legend('rhs', 'rhs new')



%%

