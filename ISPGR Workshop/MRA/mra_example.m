
%% Section 1. Simulating Data and Running the Analysis
% first, we create some simple time series that exhibit long range
% correlations. There are many methods to accomplish this, but we will use
% the method from Davies and Harte. See the fgn_sim function for details.
H = 0.9;
n = 2^13;
order = 1;
mins = 16;
maxs = floor(n/4);
scale_ratio = 2;
overlap = false;
x = fgn_sim(n, H)';
y = x + randn(n,1);


[s, fsx, fsy, f2xy, betaxy] = mra(x, y, order, mins, maxs, scale_ratio, ...
    overlap);

%% Section 2. Plotting data and performing statistial analysis of beta(s)
figure(1)
% DFA plot of X
subplot(2,2,1)
plot(logn(s, scale_ratio), logn(fsx, scale_ratio), '.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2x')
title('DFA of X Variable');
set(gca, 'FontSize', 18);

% DFA plot of Y
subplot(2,2,2)
plot(logn(s, scale_ratio), logn(fsy, scale_ratio), '.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2y')
title('DFA of Y Variable');
set(gca, 'FontSize', 18);

% DFA plot of COVxy
subplot(2,2,3)
plot(logn(s, scale_ratio), logn(f2xy, scale_ratio),'.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2xy')
title('DFA of Covariance');
set(gca, 'FontSize', 18);

% MRA plot of xy
subplot(2,2,4)
plot(s, betaxy, '.', 'MarkerSize', 20)
xlabel('s')
ylabel('beta(s)')
title('MRA of x and y');
hold on;
set(gca, 'FontSize', 18);


% asymptotic theory of multiscale regression analysis is not currently
% developed; however, we can use something called surrogation analysis in
% order construct confidence. In particular, we will make use of the
% iterated amplitude adjusted Fourier transform in order to construct
% shuffled copies of the original time series. This form of surrogate is
% appropriate for both multifractal analysis and fractal regression because
% it maintains the same power spectrum and autocorrelation as the original
% time series, while shuffling the Fourier phases. 

% This is relevant for the multifractal analysis because it
% shuffles the Fourier phases, effectively breaking any nonlinear 
% connections across time scales. This is relevant for fractal regression
% because local nonstationarities can lead to spurious correlations.
% Surrogation analysis allows us to generate. Hence, we can construct
% empirical confidence intervals for hypothesis testing while controlling
% for spurious relationships. We demonstrate the method using the simulated
% example above

% Our interest in constructing 95% confidence intervals. Hence, we will
% need at least 40 surrogates (e.g., 1/40 = 0.025, 39/40 = 0.975). 
n_surrogates = 40;

% construct n_surrogates total surrogates for each input signal
surrogate_x = iaaft(x, n_surrogates);
surrogate_y = iaaft(y, n_surrogates);

% allocate memory to store the surrogate scale-wise regression coefficients
beta_surrogates = zeros(length(s), n_surrogates);

% perform MRA on each pair of the surrogate time series. Note each
% surrogate pair corresponds to the original x and y series constructed
% above.
for i = 1:n_surrogates
    [~, ~ , ~, ~, beta_surrogates(:,i)] = mra(surrogate_x(:,i), surrogate_y(:,i),... 
    order, mins, maxs, scale_ratio, overlap);
end

% construct 95% confidence intervals as the 2.5% and 97.5% quantiles at
% each scale, s. Hence we will ultimately have length(s) confidence intervals in all, one
% for each of the original beta(s) estimated in Section !
surrogate_cis = zeros(length(s),2);
for i = 1:length(s)
   surrogate_cis(i,:) = quantile(beta_surrogates(i,:), [0.025 0.975]); 
end

% add confidence intervals to the graph of beta(s). This permits us to test
% the statistical significance at each scale of beta(s) at each scale, s.
% the graph is a bit jumpy; however, we can create loess smoothed versions
% to address that issue. If the beta(s) curve is also noisy, a loess curve
% can be used to understand the general trend. Uncomment the follow lines
% to plot CIs without smoothing. 

% subplot(2,2,4)
% plot(s, surrogate_cis(:,1), '-k');
% plot(s, surrogate_cis(:,2), '-k');

% for this example, the key things to note are that he confidence intervals
% (1) contain zero everywhere, suggesting that we have successfully
% elimiated any spurious correlation between surrogate time series, and (2)
% none of the beta(s) from the intact series fall within the confidence
% intervals of the surrogates. Thus, we can reject the null hypothesis that
% that the observed multiscale relationships result from spurious
% correlation produced from autocorrelation present in each time series.

subplot(2,2,4)
smoothed_lower_ci = fLOESS([s, surrogate_cis(:,1)], 4);
smoothed_upper_ci = fLOESS([s, surrogate_cis(:,2)], 4);
plot(s, smoothed_lower_ci, '--k');
plot(s, smoothed_upper_ci, '--k');

% Note that this is not the only pattern we can observe. In the next
% example, we examine the correlation between COP velocity in the
% anterioposterior and mediolateral directions. The idea here is that
% control processes responsible for maintaining upright posture should be
% related. 

%% Section 3. Application to postural data
copv = readtable('s007_cop.csv');
x = copv.vx;
y = copv.vy;
[s, fsx, fsy, f2xy, betaxy] = mra(x, y, order, mins, maxs, scale_ratio, ...
    overlap);
figure(2)
% DFA plot of X
subplot(2,2,1)
plot(logn(s, scale_ratio), logn(fsx, scale_ratio), '.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2x')
title('DFA of X Variable');
set(gca, 'FontSize', 18);

% DFA plot of Y
subplot(2,2,2)
plot(logn(s, scale_ratio), logn(fsy, scale_ratio), '.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2y')
title('DFA of Y Variable');
set(gca, 'FontSize', 18);

% DFA plot of COVxy
subplot(2,2,3)
plot(logn(s, scale_ratio), logn(f2xy, scale_ratio),'.', 'MarkerSize', 20)
xlabel('logs')
ylabel('logF^2xy')
title('DFA of Covariance');
set(gca, 'FontSize', 18);

% MRA plot of xy
subplot(2,2,4)
plot(s, betaxy, '.', 'MarkerSize', 20)
xlabel('s')
ylabel('beta(s)')
title('MRA of x and y');
hold on;
set(gca, 'FontSize', 18);

% Our interest in constructing 95% confidence intervals. Hence, we will
% need at least 40 surrogates (e.g., 1/40 = 0.025, 39/40 = 0.975). 
n_surrogates = 40;

% construct n_surrogates total surrogates for each input signal
surrogate_x = iaaft(x, n_surrogates);
surrogate_y = iaaft(y, n_surrogates);

% allocate memory to store the surrogate scale-wise regression coefficients
beta_surrogates = zeros(length(s), n_surrogates);

% perform MRA on each pair of the surrogate time series. Note each
% surrogate pair corresponds to the original x and y series constructed
% above.
for i = 1:n_surrogates
    [~, ~ , ~, ~, beta_surrogates(:,i)] = mra(surrogate_x(:,i), surrogate_y(:,i),... 
    order, mins, maxs, scale_ratio, overlap);
end

% construct 95% confidence intervals as the 2.5% and 97.5% quantiles at
% each scale, s. Hence we will length(s) confidence intervals in all, one
% for each of the original beta(s) estimated in Section !
surrogate_cis = zeros(length(s),2);
for i = 1:length(s)
   surrogate_cis(i,:) = quantile(beta_surrogates(i,:), [0.025 0.975]); 
end

% add confidence intervals to the graph of beta(s). This permits us to test
% the statistical significance at each scale of beta(s) at each scale, s.
% the graph is a bit jumpy; however, we can create loess smoothed versions
% to address that issue. If the beta(s) curve is also noisy, a loess curve
% can be used to understand the general trend. Uncomment the follow lines
% to plot CIs without smoothing. 

% subplot(2,2,4)
% plot(s, surrogate_cis(:,1), '-k');
% plot(s, surrogate_cis(:,2), '-k');

% for this example, the key things to note are that the confidence
% intervals are (1) not quite symmetric. This is common with empirical CIs
% and not a cause for concern. (2) Unlike the example above, it appears
% that correlations in directional velocities are not constant across
% scale. Instead, coordination seems to be confined to a region
% corresponding to relatively brief time scales. See Likens et al. (2019)
% for additional examples of that distinguish posture on teh basis of
% vision (e.g., eyes closed vs eyes open).

subplot(2,2,4)
smoothed_lower_ci = fLOESS([s, surrogate_cis(:,1)], 4);
smoothed_upper_ci = fLOESS([s, surrogate_cis(:,2)], 4);
plot(s, smoothed_lower_ci, '--k', 'LineWidth',2);
plot(s, smoothed_upper_ci, '--k', 'LineWidth',2);
