function trajInfo = generate_data_first_order_highdim(sysInfo, K, traj_input)
% T: final time
% dt: time step
% gamma, sigma: parameters
% K_handle: function handle for K(t,s), e.g., @(t,s) exp(-(t-s))

d = sysInfo.d;
N = sysInfo.N;
dt = sysInfo.dt;
gamma = sysInfo.gamma;
sigma = sysInfo.sigma;
tgrid = sysInfo.t_grid;

v = zeros(d, N);

% Initial conditions
v(:,1)      = traj_input.v0;
eta         = traj_input.noise;

for n = 2:N
    t = tgrid(n);
    memory_sum = 0;

    % Memory integral approximation via trapezoidal rule
    for m = 1:n-1
        s = tgrid(m);
        K_val = K(t, s);
        if m == 1 || m == n-1
            weight = 0.5;
        else
            weight = 1.0;
        end
        memory_sum = memory_sum + weight * K_val * v(:, m);
    end
    memory_term = dt * memory_sum;

    % Euler-Maruyama step
    v(:, n) = v(:, n-1) + (-gamma * v(:, n-1) - memory_term + eta(:, n)) * dt;
end


% for n = 2:N
%     t = tgrid(n);
%     dt = sysInfo.dt;
% 
%     % RK4 steps
%     v_prev = v(:, n-1);
%     M1 = -gamma * v_prev - memory_at(n, tgrid, K, dt, v);
%     v_temp1 = v_prev + 0.5 * dt * M1;
% 
%     M2 = -gamma * v_temp1 - memory_at(n, tgrid, K, dt, [v(:,1:n-1), v_temp1]);
%     v_temp2 = v_prev + 0.5 * dt * M2;
% 
%     M3 = -gamma * v_temp2 - memory_at(n, tgrid, K, dt, [v(:,1:n-1), v_temp2]);
%     v_temp3 = v_prev + dt * M3;
% 
%     M4 = -gamma * v_temp3 - memory_at(n, tgrid, K, dt, [v(:,1:n-1), v_temp3]);
% 
%     drift = (1/6) * (M1 + 2*M2 + 2*M3 + M4);
% 
%     % Euler update with noise
%     v(:, n) = v_prev + dt * drift + sigma * eta(:, n);
% end



trajInfo.v = v;

end


% % Define memory term function at current t and input velocity history
% function mem = memory_at(n, tgrid, K, dt, v_input)
% t = tgrid(n);
% memory_sum = 0;
% for m = 1:n-1
%     s = tgrid(m);
%     K_val = K(t, s);
%     if m == 1 || m == n-1
%         weight = 0.5;
%     else
%         weight = 1.0;
%     end
%     memory_sum = memory_sum + weight * K_val * v_input(:, m);
% end
% mem = dt * memory_sum;
% end