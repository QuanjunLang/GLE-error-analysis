function S_norm = Schur_norm(f, h)

M = 300;
N = 10000;

t = linspace(0, M, N);
dt = t(2) - t(1);

S_norm = sqrt(sum(f(t).^2 ./ h(t)).*dt);


% figure;
% plot(t, f(t).^2 ./ h(t));

end