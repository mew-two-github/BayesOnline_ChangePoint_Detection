clear; close all;
load historic.mat 
%load new.mat;
N = length(y);
%% Estimate and plot ACFs and PACFs
% Estimate ACF
[acf_y,lags,bounds] = autocorr(y,'NumLags',20);
% Plot sample ACF
figure();
subplot(1,2,1);
autocorr(y,'NumLags',20);
set(gca,'fontsize',12,'fontweight','bold');
hca = gca;
set(hca.Children(4),'LineWidth',2)
box off;
% Estimate sample PACF
pacf_y = parcorr(y,'NumLags',20);
%plot sample PACF
subplot(1,2,2);
parcorr(y,'NumLags',20);
box off;
% ACF exponentially dies, PACF goes to zero after step 1. This means that
% the model is AR(1)
%% Building the AR model
psi = y(1:N-1);
model = fitlm(psi,y(2:N),'y~x1-1');
% Residual analysis
% Residual computation
res = model.Residuals.raw;
% ACF of residuals
figure();
subplot(1,2,1); autocorr(res,'NumLags',25);
title('ACF of residuals from AR(1) model')
box off;
% PACF of residuals
subplot(1,2,2)
parcorr(res,'NumLags',25)
title('PACF of residuals from AR(1) model')
box off;
% Whiteness test
[h_model,pval_model] = lbqtest(res);
disp('Whiteness Test for Residuals results');
disp(h_model);disp(pval_model);
if h_model == 0

end
% The residuals are white
%% Load the new data
load new.mat
N = length(y);
%% Perform RLS on new data
y = y';
%res(2:N) = y(2:N)' - predict(model,y(1:N-1)');
rls_model = recursiveLS(1,model.Coefficients.Estimate);
theta = zeros(200,1);
for k = 2:200
    theta(k) = rls_model(y(k),y(k-1));
end
%% Run BoCD post 200 datapoints
% Initialise values
mu0 = 0;k0=1;lambda=50;
% N+1*N+1 matrix. Value at (r,c) is the probability value that run length
% is c-1
P_joint = zeros(N+1); P_runlength = P_joint;
P_joint(1) = 1; H = 1/lambda;
alpha = zeros(N+1-200);beta = zeros(N+1-200);
alpha(:,1) = 20; beta(:,1) = 2;
mu = zeros(N+1-200); mu(1) = mu0;
k = zeros(N+1-200); k(1) = k0;
res = zeros(N+1-200,1);
% Since MATLAB indexing starts from 1, kth row/column denotes k-1th data
% point/k-1 run length
% Implement bocd
for t = 1:N-200
    res(t) = y(t+200)-theta(t-1+200)*y(t-1+200);
    [theta(200+t),~] = rls_model(y(t+200),y(t-1+200));
    xt = res(t);
    predictive = zeros(t,1);
    % Evaluate predictive probabilities at different run lengths
    for rt = 1:t % Run length can be 0 till t
        std_dev = sqrt(beta(t,rt)*(k(rt)+1)/(alpha(t,rt)*k(rt)));
        xt_normalised = (xt-mu(rt))./std_dev;
        predictive(rt) = tpdf(xt_normalised,2*alpha(t,rt))/std_dev; 
    end
    if t ~= 1
        % Growth probabilities
        % Run length can be 1 till t
        P_joint(t,2:t+1) = (P_joint(t-1,1:t).*predictive(1:t)')*(1-H);
        % Changepoint probability
        P_joint(t,1) = sum(P_joint(t-1,1:t).*predictive(1:t)')*H;
    else
        % P(ro=0) = 1 assumed
        P_joint(2,1) = predictive(1)*(1-H);
        P_joint(1,1) = predictive(1)*H;
    end
    % Evidence
    P_evidence = sum(P_joint(t,:));
    P_joint(t,:) = P_joint(t,:)/P_evidence;
    P_runlength(t,:) = P_joint(t,:);%/P_evidence;
    % Update statements
    alpha(t+1,2:t+1) = alpha(t,1:t) + 0.5;
    beta(t+1,2:t+1) = beta(t,1:t) + k(1:t).*(xt-mu(1:t)).^2./(2*(k(1:t)+1));
    mu(2:t+1) = (k(1:t).*mu(1:t)+xt)./(k(1:t)+1);
    k(2:t+1) = k(1:t) + 1;
    %P_runlength(t,:) = P_runlength(t,:)/sum(P_runlength(t,:));
    [~, ind] = max(P_runlength(t,:));
    % Uncomment the following lines if you want to detect just the first
    % changepoint automatically.
%     if ind < t-5 % at cp, the run length will be reset
%         break
%     end
end
plot_rt_probs(P_runlength(1:t,1:t));
% By observing the probability plot we see a drastic change in run length
% at 88 where r=3 is the maximum probable scenario. So we can consider
% that indices [86,87,88] belong to the new DGP
% Also, note that plot is made after 200, so the change point is estimated
% to be at 200+86 = 286
cp  = 285;
%% Part b)
% From RLS and cp, we have
theta_initial = theta(285);
variance_initial = var(y(2:cp)-theta*y(1:cp-1));