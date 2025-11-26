function fitInfo = plot_velocity_highdim(v, sysInfo, ttl, plotON)
%PLOT_VELOCITY_HIGHDIM Plots high-dimensional velocity trajectory over time
% in linear, semilog-y, and log-log scales
%
% Inputs:
%   v [d × N]: velocity trajectory
%   sysInfo: struct with field
%       - dt: time step size
%   ttl: title for the plot

d = sysInfo.d;
t = sysInfo.t_grid;


vnorm = vecnorm(v, 2, 1);  % L2 norm at each time

if any(isnan(vnorm))
    fitInfo = [];
    return;
end
% for exponential
% idx_fit = find(t > 10 & vnorm> 1e-200);

% for powerlaw
idx_fit = find(t > 20 & vnorm> 1e-200);


t_fit = t(idx_fit);
vnorm_fit = vnorm(idx_fit);



fig = figure;



% 1. Regular scale
subplot(1,3,1);
plot(t, v', 'LineWidth', 1.0); hold on;
plot(t, vnorm, 'k--', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Velocity');
title('(Linear)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), {'||v||'}]);
grid on;

% 2. Semilogy
subplot(1,3,2);
semilogy(t, abs(v'), 'LineWidth', 1.0); hold on;
semilogy(t, vnorm, 'k--', 'LineWidth', 1.5);

% Fit log(vnorm) = a * t + b using MATLAB fit() for direct diagnostics
t_col = t_fit(:);
log_v_col = log(vnorm_fit(:));

[fitresult_semi, gof_semi] = fit(t_col, log_v_col, 'poly1');

% Extract slope and intercept
a_semilogy = fitresult_semi.p1;
b_semilogy = fitresult_semi.p2;

% Generate fitted line
vnorm_fitline_semi = exp(a_semilogy * t + b_semilogy);
semilogy(t, vnorm_fitline_semi, 'm-', 'LineWidth', 1.5);

xlabel('Time t');
ylabel('|Velocity|');
title('(Semilogy)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), ...
    {'||v||', 'Fit (semi-log)'}], 'Location', 'best');
grid on;

% 3. Log-log
subplot(1,3,3);
loglog(t(2:end), abs(v(:,2:end))', 'LineWidth', 1.0); hold on;
loglog(t(2:end), vnorm(2:end), 'k--', 'LineWidth', 1.5);

% Fit log(vnorm) = alpha * log(t) + log(C) using fit()
log_t_col = log(t_fit(:));
log_v_col = log(vnorm_fit(:));

[fitresult_loglog, gof_loglog] = fit(log_t_col, log_v_col, 'poly1');

% Extract slope and intercept
alpha_loglog = fitresult_loglog.p1;     % slope
logC_loglog = fitresult_loglog.p2;      % log-intercept
C_loglog = exp(logC_loglog);

% Fitted line: vnorm ≈ C * t^alpha
vnorm_fitline_loglog = C_loglog * t.^alpha_loglog;
loglog(t, vnorm_fitline_loglog, 'm-', 'LineWidth', 1.5);

xlabel('Time t');
ylabel('|Velocity|');
title('(Log-Log)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), ...
    {'||v||', 'Fit (log-log)'}], 'Location', 'best');
grid on;

sgtitle(ttl);


% Print diagnostics
fprintf('--- Semilogy Fit Diagnostics ---\n');
fprintf('Fitted model: log(vnorm) = %.4f * t + %.4f\n', a_semilogy, b_semilogy);
fprintf('R-squared (R^2): %.4f\n', gof_semi.rsquare);
fprintf('RMSE: %.4f\n\n', gof_semi.rmse);


% Print diagnostics
fprintf('--- Log-Log Fit Diagnostics ---\n');
fprintf('Fitted model: log(vnorm) = %.4f * log(t) + %.4f\n', alpha_loglog, logC_loglog);
fprintf('R-squared (R^2): %.4f\n', gof_loglog.rsquare);
fprintf('RMSE: %.4f\n\n', gof_loglog.rmse);

if ~ plotON
    close(fig);
end

fitInfo.a_semilogy = a_semilogy;
fitInfo.b_semilogy = b_semilogy;
fitInfo.gof_semi = gof_semi;

fitInfo.alpha_loglog = alpha_loglog;
fitInfo.logC_loglog = logC_loglog;
fitInfo.gof_loglog = gof_loglog;
end