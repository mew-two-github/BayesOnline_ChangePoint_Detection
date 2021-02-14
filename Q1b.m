clear;close all;
%% Load data
load NMRlogWell.mat;
% Plot data
plot(y);
N = length(y);
%% Initialise values
mu = 1.15;k=0.01;lambda=250;
% N+1*N+1 matrix. Value at (r,c) is the probability value that run length
% is c-1
P_joint = zeros(N+1);
P_runlength = P_joint;
P_joint(1) = 1;
H = 1/lambda;
alpha = zeros(N+1);beta = zeros(N+1);
alpha(:,1) = 20; beta(:,1) = 2;
predictive = zeros(N+1,1);
% Since MATLAB indexing starts from 1, kth row/column denotes k-1th data
% point/k-1 run length
%% Implement bocd
for t = 1:N
    xt = y(t);
    predictive = zeros(t,1);
    % Evaluate predictive probabilities at different run lengths
    for rt = 1:t % Run length can be 0 till t
        std_dev = sqrt(beta(t,rt)*(k+1)/(alpha(t,rt)*k));
        xt_normalised = (xt-mu)/std_dev;
        predictive(rt) = tpdf(xt_normalised,2*alpha(t,rt)); 
    end
    if t ~= 1
        % Growth probabilities
        for rt = 2:t+1 % Run length can be 1 till t
            P_joint(t,rt) = (P_joint(t-1,rt-1).*predictive(rt-1))*(1-H);
        end
        % Changepoint probability
        P_joint(t,1) = sum(P_joint(t-1,1:t).*predictive(1:t)')*H;
    else
        % P(ro=0) = 1 assumed
        P_joint(2,1) = predictive(1)*(1-H);
        P_joint(1,1) = predictive(1)*H;
    end
    % Evidence
    P_evidence = sum(P_joint(t,:));
    P_runlength(t,:) = P_joint(t,:)/P_evidence;
    % Update statements
    alpha(t+1,2:N+1) = alpha(t,1:N) + 0.5;
    beta(t+1,2:N+1) = beta(t,1:N) + k*(xt-mu)^2/(2*(k+1));
    mu = (k*mu+xt)/(k+1);
    k = k + 1;
    %P_runlength(t,:) = P_runlength(t,:)/sum(P_runlength(t,:));
    if isnan(sum(P_runlength(t,:)))
        break
    end
end