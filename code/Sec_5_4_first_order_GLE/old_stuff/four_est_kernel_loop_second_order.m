function stats = four_est_kernel_loop_second_order(method, all_para, B, sysInfo, kernelInfo, potentialInfo, LyapunovInfo)

perturbInfo.method = method;
result = cell(B, 1);
for b = 1:B
    fprintf('Perturbation method: \t\t%s \nTotal length of parameter: \t%d\n', method, B)    
    fprintf('Computing for the parameter: \tNo.%d\n', b);
    perturbInfo.para = all_para(b);
    result{b} = ensemble_test_S_c2_C2_second_order(sysInfo, perturbInfo, kernelInfo, potentialInfo, LyapunovInfo);
end




%%
for b = 1:B
    stats.S_norm_delta_K(b) = result{b}.S_norm_delta_K;
    stats.S_norm_tilde_K(b) = result{b}.S_norm_tilde_K;
    stats.c2(b) = result{b}.c2;
    stats.C2(b) = result{b}.C2;
end



end