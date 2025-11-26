function out = RHS(x, v, tn, v_hist, t_grid, K, V, dt, n, gamma)
    if n <= 2
        mem = 0;
    else
        t_hist = t_grid(1:n-1);
        v_past = v_hist(1:1:n-1);
        K_vals = arrayfun(@(s) K(tn, s), t_hist);
        mem = trapz(t_hist, K_vals .* v_past);
    end
    out = (-gamma * v - V(x) - mem);
end