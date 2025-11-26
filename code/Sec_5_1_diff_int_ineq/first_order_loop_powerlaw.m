close all
clear all
clc
addpaths

rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 1;

sysInfo.dt      = 0.01;          % Time step
sysInfo.T       = 100;             % Total simulation time

sysInfo.sigma   = 0;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);


all_gamma = linspace(5, 50, 15);
all_decay = linspace(1.05, 3, 20);

all_fitInfo = cell(length(all_gamma), length(all_decay));

for i = 1:length(all_gamma)
    for j = 1:length(all_decay)

        gamma = all_gamma(i);
        decay = all_decay(j);

        fprintf('gamma = %.2f, decay = %.2f\n', gamma, decay);

        amp1 = 1;
        amp2 = 1.1;
        intercept1 = 0.1;
        intercept2 = 0.1;


        threshold = 1(decay, amp1, intercept1);
        fprintf('For powerlaw kernel, the threshold for gamma is : %.2f\n', threshold)



        sysInfo.gamma = all_gamma(i);

        kernel_type = 'powerlaw';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
        kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions
        % plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log




        %% Generate trajectory
        noise = get_noise(sysInfo);

        traj_input{1}.v0 = randn(sysInfo.d, 1);
        traj_input{1}.noise = noise;

        trajInfo{1} = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});

        % traj_input{2}.v0 = randn(sysInfo.d, 1);
        % traj_input{2}.noise = noise;
        % trajInfo{2} = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});


        %%
        fitInfo = plot_velocity_highdim(trajInfo{1}.v, sysInfo, 'traj 1', false);
        % plot_velocity_highdim(trajInfo{2}.v, sysInfo, 'traj 2');
        % fitInfo = plot_velocity_highdim(abs(trajInfo{1}.v - trajInfo{2}.v), sysInfo, 'traj difference', false);

        all_fitInfo{i, j} = fitInfo;



    end
end



%% collect fitted rate
% --- Parameters ---
default_val = 0;

% --- Collect fitted rates (alpha_loglog) ---
[nGamma, nRate] = size(all_fitInfo);
rate = zeros(nGamma, nRate);

for i = 1:nGamma
    for j = 1:nRate
        try
            rate(i, j) = all_fitInfo{i, j}.alpha_loglog;
        catch
            rate(i, j) = default_val;
        end
    end
end

% --- Compute theoretical threshold ---
threshold_N = 200;
threshold_thm = zeros(1, threshold_N);
threshold_grid = linspace(all_decay(1), all_decay(end), threshold_N);
for j = 1:threshold_N
    beta = threshold_grid(j);
    threshold_thm(j) = get_theoretical_threshold_powerlaw(beta, amp1, intercept1);
end

% --- Extract R-squared values from log-log fits ---
rsquare_values = zeros(nGamma, nRate);

for i = 1:nGamma
    for j = 1:nRate
        try
            rsquare_values(i, j) = all_fitInfo{i, j}.gof_loglog.rsquare;
        catch
            rsquare_values(i, j) = default_val;
        end
    end
end

% --- Zero out rates with poor R-squared values ---
rate(rsquare_values <= 0.9900) = 0;

%%
% --- Plot results ---
figure;

% R-squared heatmap
subplot(1, 2, 1);
imagesc(all_decay, all_gamma, rsquare_values);
colorbar;
xlabel('\beta (all\_decay)');
ylabel('\gamma (all\_gamma)');
title('R-squared from Log-Log Fit');
set(gca, 'YDir', 'normal');

% Rate heatmap with threshold overlay
subplot(1, 2, 2);
imagesc(all_decay, all_gamma, -rate./all_decay);
colorbar;
xlabel('\beta (all\_decay)');
ylabel('\gamma (all\_gamma)');
title('Rate Fit');
set(gca, 'YDir', 'normal');
hold on;
plot(threshold_grid, threshold_thm, 'k-', 'LineWidth', 2);
legend('Threshold (theoretical)', 'Location', 'best');



saveas(gcf, 'power_law_decay_ex.pdf')

%%
% === Style & size presets ===
% Choose ONE of these size presets (inches). Comment out the other.
figW = 4.2; figH = 3.8;   % two-column figure (papers, LaTeX)
% figW = 3.5; figH = 3.0; % single-column compact

fontSize = 12;            % 9–10 for papers, 12–14 for slides
lineW    = 1.0;           % overlay line width

% === Figure (tight layout) ===
fig = figure('Units','inches','Position',[0 0 figW figH], 'Color','w');
% tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% === R-squared heatmap ===
% ax1 = nexttile;
% imagesc(all_decay, all_gamma, rsquare_values);
% colormap(ax1, gray);
% set(ax1,'YDir','normal','TickDir','out','FontSize',fontSize,'LineWidth',0.75,'Layer','top');
% xlabel('\beta','Interpreter','tex');       % or 'latex' if you prefer
% % ylabel('\gamma (all\_gamma)','Interpreter','tex');
% title('R^2 from Log–Log Fit','FontWeight','normal');
% colorbar; 


% === Rate heatmap + threshold overlay ===
ax2 = nexttile;
imagesc(all_decay, all_gamma, -rate./all_decay);

vals = -rate./all_decay;
vals_stretched = sign(vals - 1).*abs(vals - 1).^1.05 + 1;   % cubic stretch around 1
imagesc(all_decay, all_gamma, vals_stretched);
colormap(ax2, gray);

set(ax2,'YDir','normal','TickDir','out','FontSize',fontSize,'LineWidth',0.75,'Layer','top');
xlabel('$\beta$','Interpreter','latex');
ylabel('$a$', 'Interpreter','latex')
title('Power-law decay rate ratio $r/\beta$','FontWeight','normal','Interpreter','latex');
hold on;
plot(threshold_grid, threshold_thm, 'k-', 'LineWidth', lineW);
leg = legend('Theoretical threshold $f(\beta)$','Location','best','Interpreter','latex'); leg.Box = 'off';
colorbar; 




exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/PowerLawDecayEx.pdf', 'ContentType','vector');

