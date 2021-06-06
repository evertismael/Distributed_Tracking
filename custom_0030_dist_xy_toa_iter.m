clear; clc; close all;
addpath('Classes')  
addpath('Helpers')  
addpath('Targets')  
% -------------------------------------------------------------------------
% 1.- Define target trajetories
% -------------------------------------------------------------------------
target = car_up_right();

% -------------------------------------------------------------------------
% 2.- Generate samples of the trayectory:
% -------------------------------------------------------------------------
dt = 0.1;%1e-1;
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_iter = 1;

% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -7;

% define BSs, FC and UKF_tracker.
noise_type = 'SNR_20m'; %  SNR_center / same / SNR_20m
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();% just for reference.
fc_dist = FC_dist();% distributed ToA estimation.

% history variables:
deltas_mean_hist = zeros(scene.N_bs,N_t);
deltas_var_hist = zeros(scene.N_bs,N_t);

deltas_iter_mean_hist = zeros(scene.N_bs,N_t);
deltas_iter_var_hist = zeros(scene.N_bs,N_t);

xy_toa_hist = zeros(2,N_t);
xy_dist_toa_hist = zeros(2,N_t,scene.N_bs);
xy_dist_iter_hist = zeros(2,N_t,scene.N_bs);

% first sample initiates the tracker:
for t_idx = 1:N_t
    % read sensor
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var, act_bss] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
    [xy_dist_toa, varxy_dist_toa, bss_dist_toa_idx] = fc_dist.multilateration_toa(deltas_mean, deltas_var, act_bss);
    
    % iterations: each BS computes prior->reestimate toa-> multilateration:
    xy_dist_iter = xy_dist_toa;
    varxy_dist_iter = varxy_dist_toa;
    for iter_idx = 1:N_iter
        % prior -> refineToA -> multilateration -> saveHist:
        [prior_dist_mean, prior_dist_var, bs_from, bs_to] = fc_dist.prior_toa(xy_dist_iter, varxy_dist_iter,bss_dist_toa_idx,act_bss);
        [~, deltas_iter_mean, deltas_iter_var] = bss.refine_dist_toa(prior_dist_mean, prior_dist_var, bs_to, act_bss);
        [xy_dist_iter, varxy_dist_iter, bss_dist_toa_idx] = fc_dist.multilateration_toa(deltas_iter_mean, deltas_iter_var, act_bss);
    end
    
    % select which deltas are z:
    if N_iter > 0
        z_mean = deltas_iter_mean;
        z_var = deltas_iter_var;
    else
        z_mean = deltas_mean;
        z_var = deltas_var;
    end
    
    
    % ______________________________________________
    deltas_mean_hist(act_bss,t_idx) = deltas_mean;
    deltas_var_hist(act_bss,t_idx) = deltas_var;
    
    deltas_iter_mean_hist(act_bss,t_idx) = deltas_iter_mean;
    deltas_iter_var_hist(act_bss,t_idx) = deltas_iter_var;
    
    xy_toa_hist(:,t_idx) = xy_toa;
    xy_dist_toa_hist(:,t_idx,bss_dist_toa_idx) = xy_dist_toa;
    xy_dist_iter_hist(:,t_idx,bss_dist_toa_idx) = xy_dist_iter;
end
% -------------------------------------------------------------------------
% 4.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist,[0 90 0 90]);

fig2 = figure('Position',[2494 551 560 420]);
show_bss_xy_dist_toa(fig2,target,xy_dist_toa_hist,[0 90 0 90]);

fig3 = figure('Position',[3122 523 560 420]);
show_measurements(fig3,target.t_vect, deltas_mean_hist, deltas_var_hist);

%fig4 = figure('Position',[1932 40 560 420]);


fig5 = figure('Position',[2485 32 560 420]);
show_bss_xy_dist_toa(fig5,target,xy_dist_iter_hist,[0 90 0 90]);

fig6 = figure('Position',[3096 32 560 420]);
show_measurements(fig6,target.t_vect, deltas_iter_mean_hist, deltas_iter_var_hist);