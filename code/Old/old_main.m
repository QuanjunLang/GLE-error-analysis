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
sysInfo.gamma   = 20;         % Damping coefficient
sysInfo.dt      = 0.01;          % Time step
sysInfo.T       = 10;             % Total simulation time

sysInfo.sigma   = 1e-3;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);


%% Generate two kernels to compare
kernel_type = 'exp';                    % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
kernelInfo = kernel_info(kernel_type, sysInfo);          % Generate the kernel functions
plot_kernels(kernelInfo);               % Plot the kernel functions in three panels: linear, semilog-y, and log-log


%% Potential settings
potential_type = 'quadratic';
u = 2;

potentialInfo   = potential_info(potential_type, sysInfo, 'x_min', [1, 1]', 'u', u);
plot_potential_2d(potentialInfo, [-3, 5], [-2, 6], 200);


%% Generate trajectory
trajInfo = generate_data_second_order_highdim(sysInfo, kernelInfo.K1, potentialInfo);

plot_trajectory_highdim(trajInfo, sysInfo)
%%
% Initialize solution arrays
x1 = zeros(1, N); v1 = zeros(1, N);  % Trajectory for Kernel 1
x2 = zeros(1, N); v2 = zeros(1, N);  % Trajectory for Kernel 2

x0 = 1;
v0 = 2;

x1(1) = x0;
x2(1) = x0;

v1(1) = v0;
v2(1) = v0;
% Euler-Maruyama Integration with Memory Kernels (Synchronized Coupling)

E_difference_int = zeros(1, N);
for n = 2:N
    % Compute memory term as convolution
    memory_term1 = sum(K1(1:n-1) .* v1(n-1:-1:1) * dt);
    memory_term2 = sum(K2(1:n-1) .* v2(n-1:-1:1) * dt);
    
    % Compute velocity updates with external force (-dU/dx = -k_potential*x)
    v1(n) = v1(n-1) + dt/mass * (-gamma * v1(n-1) - memory_term1 - m*(x1(n-1)-c) + eta(n));
    v2(n) = v2(n-1) + dt/mass * (-gamma * v2(n-1) - memory_term2 - m*(x2(n-1)-c) + eta(n));
    
    % Update position
    x1(n) = x1(n-1) + v1(n) * dt;
    x2(n) = x2(n-1) + v2(n) * dt;

    E_difference_int(n) = sum(abs(K1(1:n-1) .* v1(n-1:-1:1) - K2(1:n-1) .* v2(n-1:-1:1)) * dt);
end

% Plot Memory Kernels
figure;
subplot(3,2,1);
plot(t, K1, 'b', 'LineWidth', 1.5, 'DisplayName', 'Memory Kernel 1 (\tau_1 = 1)');
hold on;
plot(t, K2, 'r', 'LineWidth', 1.5, 'DisplayName', 'Memory Kernel 2 (\tau_2 = 2)');
xlabel('Time'); ylabel('Memory Kernel K(t)');
title('Comparison of Memory Kernels');
legend();
grid on;

% Plot Results: Trajectory Comparison
subplot(3,2,3);
plot(t, x1, 'b', 'DisplayName', 'Trajectory (Kernel 1)'); hold on;
plot(t, x2, 'r', 'DisplayName', ['Traje' ...
    'ctory (Kernel 2)']);
xlabel('Time'); ylabel('Position x(t)');
title('Comparison of Trajectories with Synchronized Noise');
legend();
grid on;

subplot(3,2,5);
plot(t, v1, 'b', 'DisplayName', 'Velocity (Kernel 1)'); hold on;
plot(t, v2, 'r', 'DisplayName', 'Velocity (Kernel 2)');
xlabel('Time'); ylabel('Velocity v(t)');
title('Velocity Evolution with Different Memory Effects');
legend();
grid on;


subplot(3,2,2);
plot(t, abs(K1-K2), 'b', 'LineWidth', 1.5, 'DisplayName', 'Memory Kernel 1 (\tau_1 = 1)');
xlabel('Time'); ylabel('Memory Kernel K(t)');
title('Comparison of Memory Kernels');
legend();
grid on;

% Plot Results: Trajectory Comparison
subplot(3,2,4);
plot(t, abs(x1-x2), 'b', 'DisplayName', 'Trajectory difference'); hold on;
xlabel('Time'); ylabel('Position x(t)');
title('Comparison of Trajectories with Synchronized Noise');
legend();
grid on;


subplot(3,2,6);
plot(t, abs(v1-v2), 'b', 'DisplayName', 'Velocity difference');
xlabel('Time'); ylabel('Velocity v(t)');
title('Velocity Evolution with Different Memory Effects');
legend();
grid on;


%%
% Compute differences in position and velocity
Dx = x1 - x2;   % Difference in position
Dv = v1 - v2;   % Difference in velocity

% Compute energy-like function (relative energy between two trajectories)
E = 0.5 * Dv.^2 + 0.5 * m * Dx.^2; % Energy-like quantity

% Compute energy dissipation rate (gradient of energy)
dE_dt = gradient(E, dt); % Numerical time derivative of energy

% Compute theoretical energy dissipation term
E_dissipation = -gamma * abs(Dv).^2 + E_difference_int .* abs(Dv); % Energy balance term

% Compute cumulative integral of the dissipation term for verification
E_integral = cumsum(E_dissipation) * dt;

% Create a figure to visualize differences
figure;

% --- Plot Position and Velocity Differences ---
subplot(3, 1, 1);
plot(t, Dx, 'b', 'LineWidth', 1.5); % Blue line for position difference
ylabel('Dx = x_1 - x_2');
hold on;
plot(t, Dv, 'r', 'LineWidth', 1.5, 'LineStyle', '--'); % Red dashed line for velocity difference
ylabel('Dv = v_1 - v_2');
xlabel('Time (t)');
title('Position and Velocity Differences Between Two Trajectories');
legend('Dx(t) (Position Difference)', 'Dv(t) (Velocity Difference)');
grid on;

% --- Plot Energy Dissipation Rate (Derivative of Energy) ---
subplot(3, 1, 2);
hold on;
plot(t, dE_dt, 'b', 'LineWidth', 1.5, 'LineStyle', '-.'); % Blue dash-dot for numerical energy derivative
plot(t, E_dissipation, 'r', 'LineWidth', 1.5, 'LineStyle', '--'); % Red dashed line for theoretical dissipation
xlabel('Time (t)');
ylabel('Energy Rate');
title('Energy Dissipation Rate');
legend('dE/dt (Numerical Derivative)', '-\gamma |Dv|^2 + E_{diff} |Dv| (Theoretical)');
grid on;

% --- Plot Energy Evolution and its Theoretical Estimate ---
subplot(3, 1, 3);
hold on;
plot(t, E, 'k', 'LineWidth', 1.5, 'LineStyle', '-'); % Black solid line for energy function
plot(t, E_integral, 'r', 'LineWidth', 1.5, 'LineStyle', '--'); % Red dashed for cumulative dissipation integral
xlabel('Time (t)');
ylabel('Energy E(t)');
title('Energy Evolution');
legend('E(t) = 1/2 Dv^2 + (m/2) Dx^2', 'Cumulative Energy Dissipation');
grid on;

% Overall figure adjustments
sgtitle('Comparison of Two Trajectories with Different Memory Kernels');
