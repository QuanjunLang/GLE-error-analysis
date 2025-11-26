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

%% New figure Oct 3

ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 8 2.5]); % wide figure, short height
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

for i = 1:4
    nexttile;

    % S_norms = stats{i}.S_norm_delta_K.^2 .* stats{i}.S_norm_tilde_K.^2;
    S_norms = stats{i}.S_norm_delta_K.^2;
    
    plot(all_para{i}, stats{i}.C2./S_norms);




    % scatter(S_norms, stats{i}.C2./S_norms, 20, 'filled', 'MarkerFaceAlpha',0.6);
    % xlim([0, max(S_norms)]);
    % ylim([0, max(stats{i}.C2./S_norms)]);
    % axis square; % <-- make subplot square
    
    title(ttls{i},  'FontSize',12);
    xlabel('alpha', 'Interpreter','latex');
    

    ylabel('$C_2$', 'Interpreter','latex');

    
    grid on;
    box on;
end

set(gcf, 'Color', 'w');


%%

ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 8 2.5]); % wide figure, short height
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

for i = 1:4
    nexttile;

    % S_norms = stats{i}.S_norm_delta_K.^2 .* stats{i}.S_norm_tilde_K.^2;
    S_norms = stats{i}.S_norm_delta_K.^2;
    
    % plot(all_para{i}, stats{i}.C2./S_norms);




    scatter(S_norms, stats{i}.C2, 20, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0, max(S_norms)]);
    ylim([0, max(stats{i}.C2)]);
    axis square; % <-- make subplot square
    
    title(ttls{i},  'FontSize',12);
    xlabel('$\|K - \tilde K\|^2_{S_h}$', 'Interpreter','latex');
    

    ylabel('$C_2$', 'Interpreter','latex');

    
    grid on;
    box on;
end

set(gcf, 'Color', 'w');



%%
figure('Units','inches','Position',[1 1 8 2.5]); % wide figure, short height
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

for i = 1:4
    nexttile; hold on; grid on;
    S = stats{i}.S_norm_tilde_K.^2 .* stats{i}.S_norm_delta_K.^2;
    
    plot(stats{i}.S_norm_delta_K.^2, 'LineWidth',2.0, ...
        'DisplayName','$\|\Delta K\|^2$');
    plot(stats{i}.S_norm_tilde_K.^2, 'LineWidth',2.0, ...
        'DisplayName','$\|\tilde K\|^2$');
    plot(S, 'LineWidth',2.0, ...
        'DisplayName','$\|\tilde K\|^2 \|\Delta K\|^2$');
    plot(stats{i}.C2, 'LineWidth',2.0, ...
        'DisplayName','$C_2$');
    plot(stats{i}.C2 ./ S, 'LineWidth',2.0, ...
        'DisplayName','$C_2 / (\|\tilde K\|^2 \|\Delta K\|^2)$');
    plot(stats{i}.C2 ./ stats{i}.S_norm_tilde_K.^2, 'LineWidth',2.0, ...
        'DisplayName','$C_2 / \|\tilde K\|^2$');
    
    axis tight
    set(gca,'FontSize',14)  % larger axis labels/ticks
    
    if i == 1
        legend('Interpreter','latex','Location','best','FontSize',12);
    end
end

%% Old figure Oct 2 (Linear)
ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 8 2.5]); % wide figure, short height
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

for i = 1:4
    nexttile;
    scatter(stats{i}.S_norm.^2, stats{i}.C2, 20, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0, max(stats{i}.S_norm.^2)]);
    ylim([0, max(stats{i}.C2)]);
    axis square; % <-- make subplot square
    
    title(ttls{i},  'FontSize',12);
    xlabel('$\|K - \tilde K\|^2_{S_h}$', 'Interpreter','latex');
    
    if i == 1
        ylabel('$C_2$', 'Interpreter','latex');
    else
        ylabel('');
        set(gca,'YTickLabel',[]);
    end
    
    grid on;
    box on;
end

set(gcf, 'Color', 'w');

%%

% save('kernel_perturbation_results.mat');
% exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/Kernel_perturbation_first_order.pdf', 'ContentType','vector');


%%




% f2_info.method = 'shift';
% all_para = linspace(1, 3, B);
% stats{1} = four_est_kernel_loop(method, all_para);
% 
% 
% f2_info.method = 'rate_plus';
% all_para = linspace(0.1, 1, B);
% stats{2} = four_est_kernel_loop(method, all_para);
% 
% 
% f2_info.method = 'cut';
% all_para = linspace(0.1, 3, B);
% stats{3} = four_est_kernel_loop(method, all_para);
% 
% 
% f2_info.method = 'osc';
% all_para = linspace(0.1, 4, B);
% stats{4} = four_est_kernel_loop(method, all_para);
% 
% 
% 
% 




% 
% 
% 
% 
% for b = 1:B
%     b
% 
%     f2_info.para = all_para(b);
%     result{b} = ensemble_test_S_c2_C2(sysInfo, f2_info);
% end
% 
% 
% 
% 
% %%
% for b = 1:B
%     S_norm(b) = result{b}.S_norm;
%     c2(b) = result{b}.c2;
%     C2(b) = result{b}.C2;
% end
% 
% 
% 
% %%
% figure;hold on;
% % subplot(121);
% % scatter(S_norm.^2, c2/(sysInfo.sigma^2));
% % xlim([0, max(S_norm.^2)])
% % ylim([0, max(c2/(sysInfo.sigma^2))])
% % title('c2')
% % subplot(122);
% scatter(S_norm.^2, C2);
% xlim([0, max(S_norm.^2)])
% ylim([0, max(C2)])
% title('C2')
% % legend()