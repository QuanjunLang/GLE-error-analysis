function fitInfo = fit_traj_rate(v, sysInfo, fitOptions)
%PLOT_VELOCITY_HIGHDIM Plots high-dimensional velocity trajectory over time
% in linear, semilog-y, and log-log scales
%
% Inputs:
%   v [d × N]: velocity trajectory
%   sysInfo: struct with field
%       - dt: time step size
%   ttl: title for the plot


%% Extract Parameters
d = sysInfo.d;
t = sysInfo.t_grid;
ttl = fitOptions.ttl;
plotON = fitOptions.plotON;

fit_type = fitOptions.type;
if d  == 1
    vnorm = abs(v);
else

    vnorm = vecnorm(v, 2, 1);  % L2 norm at each time
end

if any(isnan(vnorm))        % If traj blows up, return []
    fitInfo = [];
    return;
end

%% Determine the fitting range
fit_t0 = fitOptions.fit_t0;
fit_v0 = max(fitOptions.fit_v0, min(vnorm));

idx_fit = find(t > fit_t0 & vnorm> fit_v0);
t_fit = t(idx_fit);
vnorm_fit = vnorm(idx_fit);



%%

t_col       = t_fit(:);                   % Fit log(vnorm) = a * t + b using MATLAB fit() for direct diagnostics
log_v_col   = log(vnorm_fit(:));
[fitresult_semi, gof_semi] = fit(t_col, log_v_col, 'poly1');

a_semilogy = fitresult_semi.p1;     % Extract slope and intercept
b_semilogy = fitresult_semi.p2;

fitInfo.a_semilogy = a_semilogy;
fitInfo.b_semilogy = b_semilogy;
fitInfo.gof_semi = gof_semi;




log_t_col = log(t_fit(:));          % Fit log(vnorm) = alpha * log(t) + log(C) using fit()
log_v_col = log(vnorm_fit(:));
[fitresult_loglog, gof_loglog] = fit(log_t_col, log_v_col, 'poly1');

alpha_loglog = fitresult_loglog.p1;    % Extract slope and intercept
logC_loglog = fitresult_loglog.p2;
C_loglog = exp(logC_loglog);

fitInfo.alpha_loglog = alpha_loglog;
fitInfo.logC_loglog = logC_loglog;
fitInfo.gof_loglog = gof_loglog;


if fitOptions.displayOn
    switch fit_type
        case 'semilog'          % Fit exponential (semi-log fit)
            % Print diagnostics
            fprintf('--- Semilogy Fit Diagnostics ---\n');
            fprintf('Fitted model: log(vnorm) = %.4f * t + %.4f\n', a_semilogy, b_semilogy);
            fprintf('R-squared (R^2): %.4f\n', gof_semi.rsquare);
            fprintf('RMSE: %.4f\n\n', gof_semi.rmse);

        case 'loglog' % Fit powerlaw (log-log fit)
            % Print diagnostics
            fprintf('--- Log-Log Fit Diagnostics ---\n');
            fprintf('Fitted model: log(vnorm) = %.4f * log(t) + %.4f\n', alpha_loglog, logC_loglog);
            fprintf('R-squared (R^2): %.4f\n', gof_loglog.rsquare);
            fprintf('RMSE: %.4f\n\n', gof_loglog.rmse);
    end
end

%% Plot

fig = figure;
y_min = min(vnorm)*0.1;
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

vnorm_fitline_semi = exp(a_semilogy * t + b_semilogy);      % Generate fitted line
semilogy(t, vnorm_fitline_semi, 'm-', 'LineWidth', 1.5);
xline(fit_t0)
yline(fit_v0)

xlabel('Time t');
ylabel('|Velocity|');
ylim([y_min, 10])
title('(Semilogy)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), {'||v||', 'Fit (semi-log)'}], 'Location', 'best');
grid on;



% 3. Log-log
subplot(1,3,3);
loglog(t(2:end), abs(v(:,2:end))', 'LineWidth', 1.0); hold on;
loglog(t(2:end), vnorm(2:end), 'k--', 'LineWidth', 1.5);

vnorm_fitline_loglog = C_loglog * t.^alpha_loglog;      % Fitted line: vnorm ≈ C * t^alpha
loglog(t, vnorm_fitline_loglog, 'm-', 'LineWidth', 1.5);
xline(fit_t0)
yline(fit_v0)

xlabel('Time t');
ylabel('|Velocity|');
ylim([y_min, 10])
title('(Log-Log)');
legend([arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false), {'||v||', 'Fit (log-log)'}], 'Location', 'best');
grid on;

sgtitle(ttl);

if ~ plotON
    close(fig);
end

% loglog(t, fitOptions.h(t))
end