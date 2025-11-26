function fit_result = fit_new_powerlaw(sysInfo, fitInfo)

Rt  = fitInfo.Rt;
ht  = fitInfo.ht;
c2  = fitInfo.c2;
C2  = fitInfo.C2;
t   = sysInfo.t_grid;


t0 = fitInfo.t0;
v0 = fitInfo.v0;

% Fitting range (e.g. late time, non-trivial decay)
idx_fit = find(t > t0 & Rt' > v0);
t_fit = t(idx_fit);
v_fit = Rt(idx_fit);

% 
% [fit_semi_v, gof_semi_v] = fit(t_fit(:), log(v_fit(:)), 'poly1');
% v_fit_semi = exp(fit_semi_v.p1 * t + fit_semi_v.p2);



log_t_col = log(t_fit(:));          % Fit log(vnorm) = alpha * log(t) + log(C) using fit()
log_v_col = log(v_fit(:));

[fit_semi_v, gof_semi_v] = fit(log_t_col, log_v_col, 'poly1');
v_fit_semi = exp(fit_semi_v.p2) * t.^fit_semi_v.p1;


if fitInfo.plotON

    figure;
    subplot(221);
    plot(sysInfo.t_grid, Rt./(ht + c2))
    
    subplot(222);hold on;grid on;
    loglog(sysInfo.t_grid, ht, 'DisplayName','h');

    loglog(sysInfo.t_grid, C2*(ht + c2), 'DisplayName','upper bound')
    loglog(sysInfo.t_grid, Rt, 'DisplayName','ER_t');
    set(gca, 'YScale', 'log')
    legend()


    subplot(223);hold on;
    loglog(t, Rt, 'LineWidth', 1.0, 'DisplayName', 'Lyapnuov');
    loglog(t, v_fit_semi, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Fit');
    loglog(t, C2*(ht + c2), 'm-', 'LineWidth', 1.5, 'DisplayName', 'ht');
    % loglog(t, v_fit_semi, 'm-', 'LineWidth', 1.5);
    xlabel('Time t'); ylabel('|v|');
    title('Velocity (Semilogy)');
    legend('Location', 'best'); grid on;
    set(gca, 'Yscale', 'log')
    set(gca, 'Xscale', 'log')


    subplot(224);hold on;
    plot(t, Rt, 'LineWidth', 1.0, 'DisplayName', 'Lyapnuov');
    v_fit = exp(fit_semi_v.p1 * t + fit_semi_v.p2);
    plot(t, v_fit, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Fit');






end

if fitInfo.printON
    % Print diagnostics
    fprintf('Semilog fitting: log(traj) = %.4f * t + %.4f, R^2 = %.4f\n', ...
        fit_semi_v.p1, fit_semi_v.p2, gof_semi_v.rsquare);
end







fit_result.fit_semi_v = fit_semi_v;
fit_result.gof_semi_v = gof_semi_v;



end