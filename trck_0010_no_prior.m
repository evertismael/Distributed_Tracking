clear; clc; close all;
addpath('Classes')  
addpath('Helpers')  
% -------------------------------------------------------------------------
% 1.- Define target trajetories
% -------------------------------------------------------------------------
target = oval_trayectory();

% -------------------------------------------------------------------------
% 2.- Generate samples of the trayectory:
% -------------------------------------------------------------------------
dt = (1e-1);
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);


% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -5;


% define BSs, FC and UKF_tracker.
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();

alpha = 0.2;
kappa = 3-4; % 3-n 
beta = 2;
ukf_fcnp = UKF_FcNp(alpha, beta, kappa);

% history variables:
deltas_mean_hist = zeros(scene.N_bs,N_t);
deltas_var_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
xy_ukf_hist = zeros(4,N_t);
eig_P_est_hist = zeros(4,N_t);
eig_P_pred_hist = zeros(4,N_t);


% first sample initiates the tracker:
for t_idx = 1:N_t
    %target.history(:,t_idx)
    noise_type = 'SNR_center'; %  SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var,act_bss] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
    
    if t_idx == 1
        % initialization:
        ukf_fcnp = ukf_fcnp.set_x0(xy_toa);
    else
        % tracking:
        dt = target.t_vect(t_idx) - target.t_vect(t_idx-1); % in case is not continuous sampling
        ukf_fcnp = ukf_fcnp.predict(dt);
        ukf_fcnp = ukf_fcnp.correct(deltas_mean,deltas_var);
    end
    % ______________________________________________
    % save for history:
    deltas_mean_hist(:,t_idx) = deltas_mean;
    deltas_var_hist(:,t_idx) = deltas_var;
    xy_toa_hist(:,t_idx) = xy_toa;
    xy_ukf_hist(:,t_idx) = ukf_fcnp.x_est;
    eig_P_est_hist(:,t_idx) = ukf_fcnp.eig_P_est;
    eig_P_pred_hist(:,t_idx) = ukf_fcnp.eig_P_pred;

end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------

fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist,[0 50 0 50]);
fig2 = figure('Position',[2577 858 560 420]);
show_target_and_tracker(fig2, target, xy_ukf_hist,eig_P_est_hist,eig_P_pred_hist);
fig3 = figure('Position',[1946 50 560 420]);
show_measurements(fig3, target.t_vect, deltas_mean_hist, deltas_var_hist);
''
