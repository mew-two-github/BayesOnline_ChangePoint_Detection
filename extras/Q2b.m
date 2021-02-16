%% Part b)
% From RLS and cp, we have
theta_initial = theta(285);
variance_initial = var(y(2:cp)-theta_initial*y(1:cp-1));
% Post cp
data = y(286:N);
l =  length(data);
parcorr(data,'NumLags',20);
% Again AR(1) model is sufficient to explain. This is as expected because
% only the white noise term's variance changes.
model2 = fitlm(data(1:l-1),data(2:l),'y~x1-1');
res2 = model2.Residuals.Raw;
theta_new = model2.Coefficients.Estimate;
variance_later= var(res2);