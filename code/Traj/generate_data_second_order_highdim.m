function trajInfo = generate_data_second_order_highdim(sysInfo, K, potentialInfo, traj_input)
%GENERATE_DATA_SECOND_ORDER_HIGHDIM Simulates high-dim 2nd-order Langevin system with memory
%
% Equation:
%   dX_t = V_t * dt
%   dV_t = -gamma * V_t * dt
%         - u * grad_U(X_t) * dt
%         - (integral_0^t K(t, s) * V_s ds) * dt
%         + sigma * dB_t
%
% Inputs:
%   sysInfo: struct with fields
%       - d, N, dt, gamma, u, sigma
%       - x0 [d×1], v0 [d×1]
%       - eta [d×N] Brownian increments
%   kernelInfo: struct with function handle K(t_idx, s_idx) → [d×d]
%   potentialInfo: struct with grad_U(x) → [d×1]
%
% Output:
%   trajInfo: struct with fields x [d×N], v [d×N]

% Unpack system parameters
d = sysInfo.d;
N = sysInfo.N;
dt = sysInfo.dt;
gamma = sysInfo.gamma;
tgrid = sysInfo.t_grid;


u = potentialInfo.u;
U_grad = potentialInfo.U_grad;


x = zeros(d, N);
v = zeros(d, N);

% Initial conditions
x(:,1) = traj_input.x0;
v(:,1) = traj_input.v0;
eta = traj_input.noise;

% Time stepping
for n = 2:N
    memory_term = zeros(d,1);
    

    for k = 1:n-1
        Kmat = K(tgrid(n), tgrid(k));           % [d×d] matrix
        memory_term = memory_term + Kmat * v(:,k) * dt;
    end

    % Update velocity
    v(:,n) = v(:,n-1) + dt * (...
        -gamma * v(:,n-1) ...
        - u * U_grad(x(:,n-1)) ...
        - memory_term ...
        + eta(:,n));

    % Update position
    x(:,n) = x(:,n-1) + v(:,n) * dt;
end

% Output
trajInfo.x = x;
trajInfo.v = v;
end
