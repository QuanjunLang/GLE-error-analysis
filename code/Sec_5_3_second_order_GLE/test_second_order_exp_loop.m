close all
clear all
clc
addpaths
rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 1;
sysInfo.gamma   = 10;            % Damping coefficient
sysInfo.dt      = 0.02;          % Time step
sysInfo.T       = 50;            % Total simulation time
sysInfo.sigma   = 1e-4;
sysInfo.M       = 2;
sysInfo         = update_sysInfo(sysInfo);

% kernel settings
kernelInfo.type = 'exp_test';
kernelInfo.decay = 0.5;
kernelInfo.amp = -1;
kernelInfo.intercept = 0;
kernelInfo.h = @(t) exp(-(kernelInfo.decay*2-0.1)*t);
kernelInfo.h_rate = kernelInfo.decay*2-0.1;


% potential settings
potentialInfo.type          = 'quadratic_new';
potentialInfo.x_min_temp    = zeros(sysInfo.d, 1); % This is a fake x_min G where U = xRx + G. Will need to fit the true min
potentialInfo.u             = 10;    % u is just a parameter to adjust the Lyapunov function
potentialInfo.kappa         = 10;    % kappa is the minimal eigenvalue of the Hessian of the Potential

potentialInfo   = potential_info(sysInfo, potentialInfo);


% Lyapunov settings
LG = potentialInfo.LG;
u = potentialInfo.u;
kappa = potentialInfo.kappa;
gamma = sysInfo.gamma;
p1 = LG * u * gamma^(-2);
fprintf('Requirement equation (4.7), we need L_G u gamma^(-2) < 3/4, where the actual parameter is %.4f\n', p1);
p2 = kappa * u * gamma^(-2)/2;
lambda = min(p2, 1/8);
fprintf('Requirement equation (4.4), lambda < min(1/8, kappa u gamma^(-2)/2), where the actual parameter is %.6f and lambda is %.4f \n', p2, lambda);
LyapunovInfo = Lyapunov_info(sysInfo, potentialInfo, lambda);


% fit Info
fitInfo.t0 = 10;
fitInfo.v0 = sysInfo.sigma^3+eps;
fitInfo.plotON = false;
fitInfo.printON = true;


%% Generate two kernels to compare
B = 3;
all_method = {'shift', 'rate_plus', 'cut', 'osc'};
all_para = {linspace(0, 0.5, B), linspace(0, 1, B), linspace(1, 20, B), linspace(0, 1, B)};



Num_method = length(all_method);
result = cell(Num_method, B);

for i = 1:Num_method
    perturbInfo.method = all_method{i};
    fprintf('Perturbation method: \t\t%s\n', perturbInfo.method);

    for b = 1:B
        fprintf('Computing for the parameter: \tNo.%d out of %d\n', b, B);
        perturbInfo.para = all_para{i}(b);
        result{i, b} = ensemble_test_second_order(sysInfo, perturbInfo, kernelInfo, potentialInfo, LyapunovInfo, fitInfo);
    end
end

% %%
% save('exp_data_Oct17.mat');
% %%
% load('exp_data_Oct17.mat');
%%
all_C2 = zeros(Num_method, B);
all_Se = zeros(Num_method, B);

for i = 1:Num_method
    for b = 1:B
        all_C2(i, b) = result{i, b}.C2;
        all_Se(i, b) = result{i, b}.S_norm_delta_K;
    end
end




%% Plot for C2 decreasing

ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

figure('Units','inches','Position',[1 1 10 3]); % wide figure, short height
tiledlayout(1,Num_method, 'TileSpacing','compact', 'Padding','compact');

for i = 1:Num_method
    nexttile;

    S_norms = all_Se(i, :).^2;

    scatter(S_norms, all_C2(i, :), 20, 'filled', 'MarkerFaceAlpha',0.6);
    xlim([0, max(S_norms)]);
    ylim([0, max(all_C2(i, :))]);
    axis square; % <-- make subplot square

    title(ttls{i},  'FontSize',12);
    xlabel('$|\!|\!|K - \tilde{K} |\!|\!|^2_h$', 'Interpreter','latex');

    if i == 1
        ylabel('$M$', 'Interpreter','latex');
    end


    grid on;
    box on;
    ax = gca;
    ax.FontSize = 12;
    % ax.LineWidth = 1;
end


set(gcf, 'Color', 'w');
% exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/Kernel_perturbation_second_order.pdf', 'ContentType','vector');




%% PLot for trajectory fitting

% Titles for each method
ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};

% Larger, wider figure
figure('Units','inches','Position',[0.5 0.5 13 3.8]);

% Nice readable defaults (optional)
set(groot,'defaultAxesFontSize',16);
set(groot,'defaultLegendFontSize',16);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

tl = tiledlayout(1, Num_method, 'TileSpacing','compact', 'Padding','compact');

selected_b_index = [7,7,7,8];

% Handles for shared legend
hExp = []; hFit = []; hUb = [];

% Cache the common time grid and its limits
t_full = sysInfo.t_grid;
xlims  = [t_full(1), t_full(end)];

for i = 1:Num_method
    b = selected_b_index(i);

    ax = nexttile;

    fit_semi_v = result{i, b}.fit_result.fit_semi_v;
    t   = sysInfo.t_grid;                  % ensure same grid for all plots
    v_fit_semi = exp(fit_semi_v.p1 * t + fit_semi_v.p2);

    Rt = result{i, b}.Rt;
    C2 = result{i, b}.C2;
    c2 = result{i, b}.c2;
    ht = kernelInfo.h(t);

    % Strings (upper-bound rate is common across panels)
    upbd_lgd = sprintf('upper bound $h(t)+c$, decay rate: %.2f', -abs(kernelInfo.h_rate));

    hold(ax,'on'); grid(ax,'on');

    if i == 1
        % Capture handles for the shared legend (generic labels)
        hExp = semilogy(t, Rt, 'LineWidth', 1.2, 'DisplayName', 'Expectation of Lyapunov function');
        hFit = semilogy(t, v_fit_semi, 'LineWidth', 1.8, 'DisplayName', 'Fitted decay line');
        hUb  = semilogy(t, C2*(ht + c2), 'LineWidth', 1.8, 'DisplayName', upbd_lgd);
    else
        semilogy(t, Rt, 'LineWidth', 1.2);
        semilogy(t, v_fit_semi, 'LineWidth', 1.8);
        semilogy(t, C2*(ht + c2), 'LineWidth', 1.8);
        legend(ax,'off');
    end

    % Axes formatting
    box on;
    set(ax,'YScale','log');
    xlim(ax, xlims);                       % match t-range to t_grid exactly
    xlabel(ax, 'Time $t$');
    % ylabel(ax, '$|v|$');

    % Put *this panel's* fitted rate into the title
    % title(ax, sprintf('%s (fit rate: %.2f)', ttls{i}, fit_semi_v.p1));
    title(ax, sprintf('fitted decay rate: %.2f', fit_semi_v.p1));
end

% One merged legend for the whole figure
lg = legend([hExp hFit hUb], ...
    {'Expectation of Lyapunov function $\mathbf{E}[R_t]$', 'Fitted decay line', upbd_lgd});
lg.NumColumns = 3;
lg.Box = 'off';
lg.Layout.Tile = 'south';   % place under all tiles

set(gcf, 'Color', 'w');
% exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/traj_decay_second_order.pdf', 'ContentType','vector');

%%
% 
% 
% % plot for trajectory fittings
% ttls = {'Translation', 'Dilation', 'Cutoff', 'Oscillation'};
% figure('Units','inches','Position',[1 1 10 3]); % wide figure, short height
% tiledlayout(1,Num_method, 'TileSpacing','compact', 'Padding','compact');
% 
% 
% selected_b_index = [2,2,2,3];
% for i = 1:Num_method
%     b = selected_b_index(i);
% 
%     nexttile;
%     fit_semi_v = result{i, b}.fit_result.fit_semi_v;
%     t   = sysInfo.t_grid;
%     v_fit_semi = exp(fit_semi_v.p1 * t + fit_semi_v.p2);
%     Rt = result{i, b}.Rt;
%     C2 = result{i, b}.C2;
%     c2 = result{i, b}.c2;
%     ht = kernelInfo.h(t);
% 
%     fitted_lgd = sprintf('Fitted decay rate: %.2f', fit_semi_v.p1);
%     upbd_lgd = sprintf('upper bound h(t) + c, decay rate: -%.2f', kernelInfo.h_rate);
% 
%     hold on;grid on;
%     semilogy(t, Rt, 'LineWidth', 1.0, 'DisplayName', 'Expectation of Lyapunov function');
%     semilogy(t, v_fit_semi, 'LineWidth', 1.5, 'DisplayName', fitted_lgd);
%     semilogy(sysInfo.t_grid, C2*(ht + c2), 'LineWidth', 1.5, 'DisplayName', upbd_lgd)
%     xlabel('Time t'); ylabel('|v|');
%     legend('Location', 'best', 'Interpreter','latex');
%     set(gca, 'Yscale', 'log')
%     set(gca,'FontSize', 14)
% 
% end
% 
% 
% 


%%
% 
% clc
% clear all
% load('exp_data.mat');
% 
% 
% 
% 
% 
% result_temp = result;
% 
% 
% load('exp_data_Oct16.mat');
% all_para = {linspace(0, 0.5, B), linspace(0, 1, B), linspace(1, 20, B), linspace(0, 1, B)};
% all_method = {'shift', 'rate_plus', 'cut', 'osc'};
% result(3, :) = result_temp;
% 
% save('exp_data_Oct16_new.mat');
% 
% 
