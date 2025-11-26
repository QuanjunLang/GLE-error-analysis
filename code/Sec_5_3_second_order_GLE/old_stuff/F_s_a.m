function result = F_s_a(s, a)


result = exp(a) .* igamma(s+1, a) ./ a.^(s + 1);


end

