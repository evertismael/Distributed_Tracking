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
dt = 0.01;%(1e-1);
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_iter = 1;

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
xy_iter_hist = zeros(2,N_t,N_iter);
xy_ukf_hist = zeros(4,N_t);
eig_P_est_hist = zeros(4,N_t);
eig_P_pred_hist = zeros(4,N_t);


% first sample initiates the tracker:
for t_idx = 1:N_t
    % Sumary result:
    % The higher the sampling rate the smaller Q and the smaller the
    % increment in x(k+1). So the filter assumes the object is static.
    % Since that believe is fed as prior for toa measurements, the
    % estimation of toa assumes that the object didn't move, hece the prior
    % is incorrect and therefore the final measurements. It's a sort of 
    % confirmation bias!
    
    %target.history(:,t_idx)
    noise_type = 'SNR_center'; %  SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    if t_idx == 1
        [~, deltas_mean, deltas_var, act_bss] = bss.compute_bsbnd_toa();
        [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
        % initialization:
        ukf_fcnp = ukf_fcnp.set_x0(xy_toa);
    else
        [~, deltas_mean, deltas_var] = bss.compute_bsbnd_toa();
        [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
        
        % tracking:
        dt = target.t_vect(t_idx) - target.t_vect(t_idx-1); % in case it's not continuous sampling rate
        ukf_fcnp = ukf_fcnp.predict(dt);
        
        xy_est = [ukf_fcnp.x_est(1) ukf_fcnp.x_est(3)].';
        Pxy_est = [ukf_fcnp.P_est(1,1) ukf_fcnp.P_est(1,3); ukf_fcnp.P_pred(1,3) ukf_fcnp.P_pred(3,3)];
        [prior_mean, prior_var] = fc.prior_toa(xy_est,Pxy_est, act_bss);
        [~, deltas_mean, deltas_var] = bss.refine_toa(prior_mean, prior_var, act_bss);
    
        [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
        
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
show_target_toa_meas(fig1,target,xy_toa_hist,[0,50,0,50]);
fig2 = figure('Position',[2577 858 560 420]);
show_target_and_tracker(fig2, target, xy_ukf_hist,eig_P_est_hist,eig_P_pred_hist);
fig3 = figure('Position',[1946 50 560 420]);
show_measurements(fig3, target.t_vect, deltas_mean_hist, deltas_var_hist);
'';
