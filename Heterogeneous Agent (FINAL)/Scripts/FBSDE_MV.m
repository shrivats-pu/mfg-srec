% This script solves the McKean-Vlasov FBSDE that arises from the mean
% field version of the optimal behaviour problem in single period SREC markets
% 
close all
% GLOBAL PARAMETERS
T = 1; % length of time window of SREC system
NdT = 52; % number of time-steps - this corresponds to weekly 
dt = T/NdT;
K = 2; % number of agent types
pi_k = [0.25, 0.75]; % the proportion of agents of each type
t = [0:dt:T];
P = 1; % penalty per unit of non-compliance
R = 1;
%__________________________________________________________________________
% PARAMETERS SPECIFIC TO AGENT TYPES
h = [0.2, 0.5]; % baseline of SREC production
sigma = [0.1, 0.15]; % volatility of SREC production

zeta = [1.75, 1.25];
gamma = [1.25, 1.75];
% as the drift for X_t = (a + b EYt + c Yt) dt
a = h;
b = -P ./ gamma / sum(pi_k ./ gamma);
c = P .* (1 ./ gamma + 1 ./ zeta);
% parameters for initial distribution of agents of each type
mu0_mu = [0.6, 0.2]; 
mu0_sigma = [0.1, 0.1];
% for now, chosen to be equal 
% zeta = P/(b+c); % implies zeta and gamma of 1.5 each
% gamma = -P/b;

Ndx = 401;
xmax = 1.5*R;
xmin = -0.25;
dx = (xmax-xmin)/Ndx;
x = [xmin:dx:xmax]';



Y = NaN(length(x), NdT+1, K);
p = NaN(length(x), NdT+1, K);

%
% generate initial guess of Y, Z from 
% assuming dX_t = h dt + sigma dW_t
%
% this outputs Y, Z, and mu at each t
x_r = repmat(x,1,length(x));
x_r2 = repmat(x, 1, length(t));
t_r = repmat(t, length(x), 1);
% computes Yt = E[Y_{t+1}|F_t] given the initial points of x 
% this is the expectation with respect to the randomness of the
% brownian motion (along paths of x across time, as opposed to along
% the distribution of x for a frozen point in time (across agents)
% does this in a vectorized manner for all grid-points of x
for k = 1:K
    Y(:,:,k) = normcdf((R - (x_r2 + h(k) * (T - t_r))) ./ (sigma(k) * sqrt(T - t_r))); 
end

% obtain initial guess for state distribution
% we choose this arbitrarily
for k = 1 : K
    for n = 1 : NdT+1
        sigma_eff = sqrt(mu0_sigma(k)^2 + sigma(k)^2 *(n-1)*dt);
        p(:, n, k) =  GenDiscreteDist(x, mu0_mu(k)+ h(k)*dt*(n-1), sigma_eff);
    end
    figure(1);
    plot(x,p(:,:,k));
    hold on
    getframe();
end

% need to store a function for the 'actual' (non-discretized) distribution
% of X
% plotting the initial distributions

prev_p = p;
CDF_old = cumsum(prev_p);
muIter = 30;
tol = 5e-5;
maxIter = 100;
counter = 1;
curr_diff = 1;
%while curr_diff > tol
% could make this a while loop using the above - ultimately chose not to in
% order to avoid having an additional parameter to worry about
for M = 1:muIter
    % step back in time to update Y & Z given distribution mu
    fprintf("FBSDE given mu...")
    tic;
    % at each time-step, update the estimate of Y_t^k
    for n = NdT : -1 : 1
        % uses current estimate of Yt, as well as the distribution and the
        % value of Y_{t+1} in order to refine the guess of Y_t by stepping
        % backwards.
        % this gets the value of Y_t for all x, t, and k
        Y(:,n,:) =  StepBack_given_mu(x, squeeze(p(:,n,:)), ...
            squeeze(Y(:,n,:)), squeeze(Y(:,n+1,:)),   dt, 10, a, b, c, sigma, pi_k, gamma, K);

        % the above returns all iterations of Y_t - therefore, we set Y_t
        % equal to the final iteration
       %  Y(:,n) = Yt_all(:,end); 
    end
    % done up to here?
    toc;
    fprintf("... done FBSDE given mu\n")
    
    % now update distribution given Y which results in
    %
    % generate the distribution X for each t given that
    %   dX = ( a + b E[Y_t] + c Y_t ) dt + sigma dW_t
    % and X_0 ~ p0
    %
    fprintf("update mu...")
    tic;
    p = UpdateMu(x, p(:,1,:), Y, sigma, dt, a, b, c, gamma, pi_k);
    toc;
    fprintf("... done update mu\n")
    
    pp = zeros(length(x), NdT+1);
    for k = 1:K
        pp = pp + pi_k(k) .* p(:,:,k);
    end
    figure(2);
    plot(x,pp);
    getframe();

    figure(3);
    plot(x,p(:, end, 1) - prev_p(:, end, 1))
    hold on
    getframe();
    
    figure(4);
    plot(x,p(:, end, 2) - prev_p(:, end, 2))
    hold on
    getframe();
    
    CDF_new = cumsum(p);
    prev_diff = curr_diff;
    curr_diff = (max(max(max(abs(CDF_new - CDF_old), [], 1), [], 2)));
    figure(5)
    plot(counter, curr_diff,'^r', 'MarkerFaceColor','r')
    hold on
    getframe();
    
    figure(6)
    plot(counter, curr_diff / prev_diff, '^r', 'MarkerFaceColor', 'r')
    hold on
    getframe();
    prev_p = p;
    CDF_old = cumsum(prev_p);

    fprintf("compute strategy...")
    tic;
    
    % computes price and controls from the variational analysis equations -
    % also creates useful plots
    [S, g, Gamma] =  ComputePrice_and_Strategy(x, p, Y, zeta, gamma, h, P, t, K, pi_k);
    toc;
    fprintf("... done compute strategy\n")
    counter = counter + 1;
end

w = plotSummary(x, p, g, Gamma, S, h, t, K, pi_k);

% simulating finite player game
Nplayers = 200;
S_finite = FinitePlayerGame(x, p, Y, zeta, gamma, sigma, h,R, P, t, K, pi_k, Nplayers, mu0_mu, mu0_sigma);

% simulating path for representative agents from each sub-population
% NSims = 1;
% [bank, cash, gen, trade] =SimPath(t, x, [0.6, 0.2], g, Gamma, S, zeta, gamma, sigma, h, P, R, K, Y, NSims);
