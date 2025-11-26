close all
clear all
clc
addpaths
rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 1;
sysInfo.gamma   = 3;         % Damping coefficient
sysInfo.dt      = 0.02;          % Time step
sysInfo.T       = 40;             % Total simulation time
sysInfo.sigma   = 1e-3;

sysInfo         = update_sysInfo(sysInfo);

%% Generate two kernels to compare
B = 10;

method = 'shift';
all_para{1} = linspace(0, 3, B);
stats{1} = four_est_kernel_loop(method, all_para{1}, B, sysInfo);


method = 'rate_plus';
all_para{2} = linspace(0, 3, B);
stats{2} = four_est_kernel_loop(method, all_para{2}, B, sysInfo);


method = 'cut';
all_para{3} = linspace(1, 5, B);
stats{3} = four_est_kernel_loop(method, all_para{3}, B, sysInfo);


method = 'osc';
all_para{4} = linspace(0.1, 4, B);
stats{4} = four_est_kernel_loop(method, all_para{4}, B, sysInfo);


save('power_law_data_new_sigma.mat');
%%

load('power_law_data_new_sigma.mat');

%%

ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 10 3]); % wide figure, short height
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

for i = 1:4
    nexttile;

    S_norms = stats{i}.S_norm_delta_K.^2;

    scatter(S_norms, stats{i}.C2, 20, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0, max(S_norms)]);
    ylim([0, max(stats{i}.C2)]);
    axis square; % <-- make subplot square

    title(ttls{i},  'FontSize',12);
    xlabel('$|\!|\!|K - \tilde{K} |\!|\!|^2_h$', 'Interpreter','latex');

    if i == 1
        ylabel('$C_2$', 'Interpreter','latex');
    end

    
    grid on;
    box on;
    ax = gca;
    ax.FontSize = 12;
    % ax.LineWidth = 1;
end


set(gcf, 'Color', 'w');
exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/Kernel_perturbation_first_order.pdf', 'ContentType','vector');

%%

ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 10 3]); % wide figure, short height
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

for i = 1:4
    nexttile;

    S_norms = stats{i}.S_norm_delta_K.^2;

    scatter(all_para{i}, stats{i}.C2./S_norms, 20, 'filled', 'MarkerFaceAlpha',0.6);
    % xlim([0, max(S_norms)]);
    % ylim([0, max(stats{i}.C2)]);
    axis square; % <-- make subplot square

    title(ttls{i},  'FontSize',12);
    xlabel('$|\!|\!|K - \tilde{K} |\!|\!|^2_h$', 'Interpreter','latex');

    if i == 1
        ylabel('$C_2$', 'Interpreter','latex');
    end


    grid on;
    box on;
    set(gca, 'FontSize', 12)
end