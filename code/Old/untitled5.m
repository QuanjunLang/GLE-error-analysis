close all; clear all; clc;

%% Parameters
mass = 1;
gamma = 10;
dt = 0.002;
T = 8;
N = round(T / dt);
t = (0:N-1) * dt;

%% Choose memory kernels
kernel_type = 'exp';  % options: 'exp', 'power', 'osc'
[f1, f2, K1, K2] = switch_kernel(kernel_type);

%% Choose external potential
V_type = 'quadratic';  % options: 'quadratic', 'doublewell', 'linear', 'sinusoidal'
V = get_potential(V_type);

%% Shared noise
sigma = 0;
eta = sqrt(2 * gamma * sigma^2 / dt) * randn(1, N);

%% Initial conditions
x1_0 = 1; v1_0 = 2;
x2_0 = 0; v2_0 = 3;

x1 = zeros(1, N); v1 = zeros(1, N);
x2 = zeros(1, N); v2 = zeros(1, N);
x1(1) = x1_0; v1(1) = v1_0;
x2(1) = x2_0; v2(1) = v2_0;

%% Memory integral storage
mem1_vals = zeros(1, N);
mem2_vals = zeros(1, N);

%% Time stepping with RK4
for n = 2:N
    tn = t(n);

    % Compute memory integrals
    if n > 2
        t_hist = t(1:n-1);
        v1_hist = v1(n-1:-1:1);
        v2_hist = v2(n-1:-1:1);

        K1_vals = arrayfun(@(s) K1(tn, s), t_hist);
        K2_vals = arrayfun(@(s) K2(tn, s), t_hist);

        mem1_vals(n) = trapz(t_hist, K1_vals .* v1_hist);
        mem2_vals(n) = trapz(t_hist, K2_vals .* v2_hist);
    else
        mem1_vals(n) = 0;
        mem2_vals(n) = 0;
    end

    if n == 2
        % Use Euler step
        v1(n) = v1(n-1) + dt/mass * (-gamma * v1(n-1) - V(x1(n-1)) - mem1_vals(n) + eta(n));
        v2(n) = v2(n-1) + dt/mass * (-gamma * v2(n-1) - V(x2(n-1)) - mem2_vals(n) + eta(n));
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

%% Plot memory integrals
figure;
plot(t, mem1_vals, 'b', 'LineWidth', 2, 'DisplayName', '$\\int_0^t K_1(t,s)v_1(s)ds$');
hold on;
plot(t, mem2_vals, 'r--', 'LineWidth', 2, 'DisplayName', '$\\int_0^t K_2(t,s)v_2(s)ds$');
xlabel('Time'); ylabel('Memory Integral');
title('Memory Term Comparison');
legend('Interpreter', 'latex');
grid on;
