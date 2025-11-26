function fitInfo = plot_trajectory_decay_fit_new(trajInfo, sysInfo, potentialInfo, ttl, plotON)
%PLOT_TRAJECTORY_DECAY_FIT Visualizes and fits decay of high-dimensional x and v trajectories
%
% Inputs:
%   trajInfo: struct with fields:
%       - x [d × N]: position trajectory
%       - v [d × N]: velocity trajectory
%   sysInfo: struct with fields:
%       - t_grid [1 × N]: time grid
%   ttl: string, title for the plot
%   plotON: boolean, whether to show the plot or not
%
% Output:
%   fitInfo: struct with fit results and goodness-of-fit metrics

x = trajInfo.x;
v = trajInfo.v;
[d, N] = size(x);
t = sysInfo.t_grid;

x_min = potentialInfo.U_argmin;
xnorm = vecnorm(x - x_min, 2, 1);
vnorm = vecnorm(v, 2, 1);

% Fitting range (e.g. late time, non-trivial decay)
idx_fit = find(t > 20 & xnorm > 1e-200 & vnorm > 1e-200);
t_fit = t(idx_fit);
x_fit = xnorm(idx_fit);
v_fit = vnorm(idx_fit);

% Initialize figure
fig = figure;

% ===================== Plot Position =====================
subplot(2, 3, 1);  % Linear
plot(t, x', 'LineWidth', 1.0); hold on;
plot(t, xnorm, 'k--', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|x|');
title('Position (Linear)');
legend([arrayfun(@(i) sprintf('x_%d', i), 1:d, 'UniformOutput', false), {'||x||'}], 'Location', 'best'); grid on;

subplot(2, 3, 2);  % Semilog
semilogy(t, abs(x'), 'LineWidth', 1.0); hold on;
semilogy(t, xnorm, 'k--', 'LineWidth', 1.5);
% Fit semilogy
[fit_semi_x, gof_semi_x] = fit(t_fit(:), log(x_fit(:)), 'poly1');
x_fit_semi = exp(fit_semi_x.p1 * t + fit_semi_x.p2);
semilogy(t, x_fit_semi, 'm-', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|x|');
title('Position (Semilogy)');
legend({'x_i', '||x||', 'Fit (semi-log)'}, 'Location', 'best'); grid on;

subplot(2, 3, 3);  % Loglog
loglog(t(2:end), abs(x(:,2:end))', 'LineWidth', 1.0); hold on;
loglog(t(2:end), xnorm(2:end), 'k--', 'LineWidth', 1.5);
% Fit loglog
[fit_log_x, gof_log_x] = fit(log(t_fit(:)), log(x_fit(:)), 'poly1');
x_fit_loglog = exp(fit_log_x.p2) * t.^fit_log_x.p1;
loglog(t, x_fit_loglog, 'm-', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|x|');
title('Position (Log-Log)');
legend({'x_i', '||x||', 'Fit (log-log)'}, 'Location', 'best'); grid on;

% ===================== Plot Velocity =====================
subplot(2, 3, 4);  % Linear
plot(t, v', 'LineWidth', 1.0); hold on;
plot(t, vnorm, 'k--', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|v|');
title('Velocity (Linear)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), {'||v||'}], 'Location', 'best'); grid on;

subplot(2, 3, 5);  % Semilog
semilogy(t, abs(v'), 'LineWidth', 1.0); hold on;
semilogy(t, vnorm, 'k--', 'LineWidth', 1.5);
[fit_semi_v, gof_semi_v] = fit(t_fit(:), log(v_fit(:)), 'poly1');
v_fit_semi = exp(fit_semi_v.p1 * t + fit_semi_v.p2);
semilogy(t, v_fit_semi, 'm-', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|v|');
title('Velocity (Semilogy)');
legend({'v_i', '||v||', 'Fit (semi-log)'}, 'Location', 'best'); grid on;

subplot(2, 3, 6);  % Loglog
loglog(t(2:end), abs(v(:,2:end))', 'LineWidth', 1.0); hold on;
loglog(t(2:end), vnorm(2:end), 'k--', 'LineWidth', 1.5);
[fit_log_v, gof_log_v] = fit(log(t_fit(:)), log(v_fit(:)), 'poly1');
v_fit_loglog = exp(fit_log_v.p2) * t.^fit_log_v.p1;
loglog(t, v_fit_loglog, 'm-', 'LineWidth', 1.5);
xlabel('Time t'); ylabel('|v|');
title('Velocity (Log-Log)');
legend({'v_i', '||v||', 'Fit (log-log)'}, 'Location', 'best'); grid on;

sgtitle(ttl, 'Interpreter', 'latex');

% Store fit info
fitInfo.position.semilog = struct('a', fit_semi_x.p1, 'b', fit_semi_x.p2, 'gof', gof_semi_x);
fitInfo.position.loglog  = struct('alpha', fit_log_x.p1, 'logC', fit_log_x.p2, 'gof', gof_log_x);

fitInfo.velocity.semilog = struct('a', fit_semi_v.p1, 'b', fit_semi_v.p2, 'gof', gof_semi_v);
fitInfo.velocity.loglog  = struct('alpha', fit_log_v.p1, 'logC', fit_log_v.p2, 'gof', gof_log_v);

% Print diagnostics
fprintf('\n=== Position Decay Fit ===\n');
fprintf('Semilog: log(xnorm) = %.4f * t + %.4f, R^2 = %.4f\n', ...
    fit_semi_x.p1, fit_semi_x.p2, gof_semi_x.rsquare);
fprintf('Loglog : log(xnorm) = %.4f * log(t) + %.4f, R^2 = %.4f\n', ...
    fit_log_x.p1, fit_log_x.p2, gof_log_x.rsquare);

fprintf('\n=== Velocity Decay Fit ===\n');
fprintf('Semilog: log(vnorm) = %.4f * t + %.4f, R^2 = %.4f\n', ...
    fit_semi_v.p1, fit_semi_v.p2, gof_semi_v.rsquare);
fprintf('Loglog : log(vnorm) = %.4f * log(t) + %.4f, R^2 = %.4f\n', ...
    fit_log_v.p1, fit_log_v.p2, gof_log_v.rsquare);

if ~plotON
    close(fig);
end

end
