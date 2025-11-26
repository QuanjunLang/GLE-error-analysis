function plot_kernels(kernelInfo)
    % Extract fields from kernelInfo
    f1 = kernelInfo.f1;
    f2 = kernelInfo.f2;
    label1 = kernelInfo.label1;
    label2 = kernelInfo.label2;
    title_str = kernelInfo.title_str;
    tau_vals = linspace(1e-1, 100, 1000);

    % Evaluate kernels
    K1_vals = f1(tau_vals);
    K2_vals = f2(tau_vals);
    absK1 = abs(K1_vals);
    absK2 = abs(K2_vals);

    % Create subplots
    figure;
    fontSize = 18; % Set global font size
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Regular scale plot
    ax1 = nexttile;
    plot(tau_vals, K1_vals, 'b', 'DisplayName', label1); hold on;
    plot(tau_vals, K2_vals, 'r--', 'DisplayName', label2);
    title(title_str, 'Interpreter', 'latex');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$K(t)$', 'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'Location', 'northeast');
    grid on;
    set(ax1, 'FontSize', fontSize);

    % Semilog-y plot
    ax2 = nexttile;
    semilogy(tau_vals, absK1, 'b', 'DisplayName', '$|K_1(t)|$'); hold on;
    semilogy(tau_vals, absK2, 'r--', 'DisplayName', '$|K_2(t)|$');
    title('Semilog-y: $\log |K(t)|$ vs $t$', 'Interpreter', 'latex');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\log |K(t)|$', 'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'Location', 'northeast');
    grid on;
    set(ax2, 'FontSize', fontSize);

    % Log-log plot
    ax3 = nexttile;
    loglog(tau_vals, absK1, 'b', 'DisplayName', '$|K_1(t)|$'); hold on;
    loglog(tau_vals, absK2, 'r--', 'DisplayName', '$|K_2(t)|$');
    title('Log-Log: $|K(t)|$ vs $t$', 'Interpreter', 'latex');
    xlabel('$\log t$', 'Interpreter', 'latex');
    ylabel('$\log |K(t)|$', 'Interpreter', 'latex');
    legend('Interpreter', 'latex', 'Location', 'northeast');
    grid on;
    set(ax3, 'FontSize', fontSize);

    % Final layout
    title(t, title_str, 'Interpreter', 'none', 'FontSize', 14);

    set(gcf, 'Position', [100, 100, 1400, 400]);
end
