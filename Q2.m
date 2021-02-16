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
model = fitlm(psi,y(2:N));
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
% The residuals are white
%% Load the new data
load new.mat
N = length(y);
%% Compute residuals
res = y;res(2:N) = y(2:N);