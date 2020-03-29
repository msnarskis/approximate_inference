%% Model (v2) Predictions

% stimulus parameters
eps_range = [0, 2];
k = 3;
n = 10;

% model parameters
sig_t = 1; % (std dev)
sig_n = 1;
sig_v = 0.1;
sig_sa = 100;
sig_sv = 100;

pr_R = 0.5;
pr_C = 0.5;

nsamp = 1; % trials per W vector

% noisy stimulus generation
noise_a = 1;
noise_v = 0.1;

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
W = [W; ones(1, k); -1*ones(1, k)];
w = size(W,1);

% correct answer
corr = repmat([ones(2^k-2,1);zeros(2,1)],1,n)';

% stimulus var
stim = generate_stimulus_v2(eps_range, n, k, W);
eps = linspace(eps_range(1), eps_range(2), n);

%% Run Model

% returns resp: (n: locations, w: size(W))

center = model_v2(stim, 0, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp);

center_corr = abs(corr - center);

match = model_v2(stim, 1, sig_t^2, sig_n^2, sig_v^2, sig_sa^2, sig_sv^2,...
    pr_R, pr_C, noise_a, noise_v, nsamp);

match_corr = abs(corr - match);

%% Model Analyses

% match per W
figure;
hold on;

title("Match R=1");
xlabel("stimulus location (eps)");
ylabel("prob R=1");
ylim([0 1]);

for i = 1:w-2
    plot(eps, match_corr(:,i), 'LineWidth', 4);
end

for i = w-1:w
    plot(eps, match_corr(:,i), '--', 'LineWidth', 4);
end


% center per W
figure;
hold on;

title("Center R=1");
xlabel("stimulus location (eps)");
ylabel("prob R=1");
ylim([0 1]);

for i = 1:w-2
    plot(eps, center_corr(:,i), 'LineWidth', 4);
end

for i = w-1:w
    plot(eps, center_corr(:,i), '--', 'LineWidth', 4);
end

%% 2D of above
figure
title("Match")
imagesc(eps, 1:w, match');
xlabel("eps")
ylabel("W")
colorbar;

figure
title("Center")
imagesc(eps, 1:w, center');
xlabel("eps")
ylabel("W")
colorbar;

%% Match vs Center
figure
hold on

title("Center (b) v Matched (r)");
xlabel("stimulus location (eps)");
ylabel("prob correct");
ylim([0 1]);

for i = 1:w-2
    plot(eps, match_corr(:,i), 'r', 'LineWidth', 4);
    plot(eps, center_corr(:,i), 'b', 'LineWidth', 4);
end

for i = w-1:w
    plot(eps, match_corr(:,i), 'r--', 'LineWidth', 4);
    plot(eps, center_corr(:,i), 'b--', 'LineWidth', 4);
end












