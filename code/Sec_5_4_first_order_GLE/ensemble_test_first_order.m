function result = ensemble_test_first_order(sysInfo, perturbInfo, kernelInfo, fitInfo)


%% update the kernel info and compute kernel error
kernelInfo = new_kernel_info_structure(kernelInfo, sysInfo, perturbInfo);


%% Generate Ensemble trajectories
M = sysInfo.M;
traj_1 = zeros(M, sysInfo.d, sysInfo.tn);
traj_2 = zeros(M, sysInfo.d, sysInfo.tn);

v0_sigma = 0.5;
v0_mu = 1;

for i = 1:M

    noise = get_noise(sysInfo);

    traj_input{1}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{1}.noise = noise;

    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});
    traj_1(i, :, :) = temp.v;

    % traj_input{2}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{2}.v0 = traj_input{1}.v0;
    traj_input{2}.noise = noise;
    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});
    traj_2(i, :, :) = temp.v;
end

traj_diff = squeeze(mean((traj_1 - traj_2).^2, 1));

%% Fit and find c2 C2
d = sysInfo.d;
if d  == 1
    vnorm2 = abs(traj_diff);
else
    vnorm2 = vecnorm(traj_diff, 2, 1)';  % L2 norm at each time
end


ht = kernelInfo.h(sysInfo.t_grid)';
c2 = sysInfo.sigma.^2;
C2 = max(vnorm2./(ht + c2));
%%

fitInfo.Rt      = vnorm2;
fitInfo.ht      = ht;
fitInfo.c2      = c2;
fitInfo.C2      = C2;


fitON = ~(max(vnorm2)<fitInfo.v0);


if fitON
    fit_result = fit_new_powerlaw(sysInfo, fitInfo);
    result.fit_result = fit_result;
else
    result.fit_result = [];
end


%%
result.C2 = C2;
result.c2 = c2;
result.S_norm_delta_K = kernelInfo.S_norm_delta_K;
result.S_norm_tilde_K = kernelInfo.S_norm_tilde_K;
result.Rt = vnorm2;




%%



end




