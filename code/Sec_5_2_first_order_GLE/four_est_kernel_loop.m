function stats = four_est_kernel_loop(method, all_para, B, sysInfo)

f2_info.method = method;
for b = 1:B
    b

    f2_info.para = all_para(b);
    result{b} = ensemble_test_S_c2_C2(sysInfo, f2_info);
end




%%
for b = 1:B
    stats.S_norm_delta_K(b) = result{b}.S_norm_delta_K;
    stats.S_norm_tilde_K(b) = result{b}.S_norm_tilde_K;
    stats.c2(b) = result{b}.c2;
    stats.C2(b) = result{b}.C2;
end



end