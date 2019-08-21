function [eps_a] = generate_stimulus(mu_a, sig_a, k, w, ext_samp)

%Generate stimulus: row vector of length k where eps_a = eps_v
% and ~ w * N(mu_ahat, sig_a) + (1-w) * N(-mu_ahat, sig_a)
% where mu_ahat is s.t. the mean of eps_a is mu_a
%
% condition 0  -->  p = 1,       k = 1
% condition 1  -->  p = 1,       k < 1
% condition 2  -->  p in (.5,1], k < 1 

mu_hat = mu_a / ((2*w - 1)); % *sqrt(k)); % effective mean
eps_a_nW = randn(k, size(mu_hat, 2), ext_samp)*sqrt(sig_a) + repmat(mu_hat, [k, 1, ext_samp]); % audio cue locations for each k

W = sign(binornd(1, w, size(eps_a_nW)) - 0.5); % correct response
sum(count(string(W), '-1'), 2);

eps_a = W .* eps_a_nW;
end