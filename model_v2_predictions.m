%% Model (v2) Predictions

% stimulus parameters
eps_range = [0, 1];
k = 3;
n = 5;

% model parameters
sig_t = 0.1; % (std dev)
sig_n = 0.1;
sig_v = 0.01;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector

% noisy stimulus generation
noise_a = .1;
noise_v = .1;

%% Generate Stimulus
% for R = 0
W = ones(2^k-2, k);

for w = 1:2^k-2
    c = dec2bin(w,k);
    for i=1:k
        W(w,i) = W(w,i) * sign(str2num(c(i))-0.5);
    end
end

% for R = 1
W = [W; ones((2^k-2)/2, k); -1*ones((2^k-2)/2, k)];
w = size(W,1);

% correct answer
corr = repmat([ones(2^k-2,1);zeros(2^k-2,1)],1,n)';

% stimulus var
stim = generate_stimulus_v2(eps_range, n, k, W);
eps = abs(stim(1,:,1)');

%% Run Model

% returns resp: (n: locations, w: size(W))

match = model_v2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp);

match_corr = abs(corr - match);

center = model_v2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp);

center_corr = abs(corr - center);

%% Model Analysis

figure;
hold on;

title("Prob correct wrt eps");
xlabel("stimulus location (eps)");
ylabel("prob R=1 or correct");
ylim([0 1]);

for i = 1:w
    plot(eps, match(:,i), 'LineWidth', 2);
    plot(eps, match_corr(:,i), '--', 'LineWidth', 2);
end














