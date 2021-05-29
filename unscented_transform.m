close all; clear; clc;
addpath('Helpers')  
% define gaussian:
mu = [0,0];
Sigma = [32,15;15,40];
R = mvnrnd(mu,Sigma,100000);


% define non-linear function:
f = @(x,y)[x + y, x.^2 + 0.001*y.^2];
R_img = f(R(:,1),R(:,2));

% true:
mean_true = mean(R_img,1);
var_xx = mean((R_img-mean_true).^2,1);
var_xy = mean(prod(R_img-mean_true,2),1);
var_true = [var_xx, var_xy]

% linear:
mean_lin = f(mu(1),mu(2));

% ukf:

% generate sigma points and weights:
n = 2; alpha = 0.2;
kappa = 3-n; 
beta = 2;

%sigma points
[sig_wm, sig_wc] = sigma_weights(n,alpha, beta, kappa);
sig_x = sigma_points(alpha, kappa, mu, Sigma);

% mean:
sig_z = zeros(2,2*n+1);
for sig_idx = 1:2*n+1
  sig_z(:,sig_idx) = f(sig_x(1,sig_idx), sig_x(2,sig_idx));
end
z_ut = sum(sig_wm.*sig_z,2);

% var:
Pz_ut = zeros(size(sig_z,1),size(sig_z,1));
for sig_idx = 1:2*n+1
 dz = (sig_z(:,sig_idx) - z_ut);
 Pz_ut = Pz_ut + sig_wc(:,sig_idx)*(dz*dz.');
end
Pz_ut(:).'

figure;
subplot(1,2,1);
scatter(R(:,1),R(:,2),'+','MarkerEdgeAlpha',0.1);


subplot(1,2,2); hold on;
scatter(R_img(:,1),R_img(:,2),'+','MarkerEdgeAlpha',0.1);
scatter(mean_lin(1),mean_lin(2),'h');
scatter(mean_true(1),mean_true(2),500,'.');
scatter(z_ut(1),z_ut(2),500,'*');

ylim([0,200]);
xlim([-100, 100]);