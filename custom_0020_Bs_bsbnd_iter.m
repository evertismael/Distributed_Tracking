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
dt = 0.01;%1e-1;
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_iter = 1;

% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -15;


% define BSs, FC and UKF_tracker.
noise_type = 'SNR_20m'; %  SNR_center / same / SNR_20m
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();% just for reference.

alpha = 0.2;
kappa = 3-4; % 3-n 
beta = 2;
N_bs = size(bss.bx,2);
dukfs = DUKFs(N_bs, alpha, beta, kappa);

% history variables:
deltas_mean_hist = zeros(scene.N_bs,N_t);
deltas_var_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
xy_iter_toa_hist = zeros(2,N_t);
xy_ukf_hist = zeros(4,N_t,N_bs);
% first sample initiates the tracker:
for t_idx = 1:N_t
    % read sensor
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var, act_bss] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
    
    % iterations:
    xy_iter = xy_toa;
    varxy_iter = varxy_toa;
    for iter_idx = 1:N_iter
        % prior -> refineToA -> multilateration:
        [prior_mean, prior_var] = fc.prior_toa(xy_iter, varxy_iter, act_bss);
        [~, deltas_iter_mean, deltas_iter_var] = bss.refine_toa(prior_mean, prior_var, act_bss);
        [xy_iter, varxy_iter] = fc.multilateration_toa(deltas_iter_mean, deltas_iter_var, act_bss);
    end
    
    if N_iter == 0
        z_mean = deltas_mean;
        z_var = deltas_var;
    else
        z_mean = deltas_iter_mean;
        z_var = deltas_iter_var;
    end
    
    
    % ______________________________________________
    % save for history:
    % measurements, cenralized loc, ukfilters, eig_P
    deltas_mean_hist(act_bss,t_idx) = z_mean;
    deltas_var_hist(act_bss,t_idx) = z_var;
    xy_toa_hist(:,t_idx) = xy_toa;
    xy_iter_toa_hist(:,t_idx) = xy_iter;
    
    %[x_ests,eig_P_ests, eig_P_preds] = dukfs.data_for_history(act_bss);
    %xy_ukf_hist(:,t_idx,:) = x_ests;
end
% -------------------------------------------------------------------------
% 4.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist,[0 90 0 90]);

fig2 = figure('Position',[3275 870 560 420]);
show_target_toa_meas(fig2,target,xy_iter_toa_hist,[0 90 0 90]);

fig3 = figure('Position',[2491 854 570 413]);
show_target_and_tracker_diffusion(fig3, target, xy_ukf_hist,[0 90 0 90]);
%%
% fig = figure;
% show_target(fig,target,[0 90 0 90]);
% 
% fig1 = figure('Position', [3256 346 560 420]);
% animate_target(fig1,target,[0 90 0 90]);
