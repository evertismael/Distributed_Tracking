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
SNR_db = -5;

% define BSs, FC and UKF_tracker.
noise_type = 'SNR_20m'; %  SNR_center / same / SNR_20m
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();% just for reference.
fcd = FC_dist();% distributed ToA estimation.

% history variables:
z_mean_hist = zeros(scene.N_bs,N_t);
z_var_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
xyd_toa_hist = zeros(2,N_t,scene.N_bs);

% first sample initiates the tracker:
for t_idx = 1:N_t
    % read sensor
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var, act_bss] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var, act_bss);
    [xyd_toa, varxyd_toa, bss_idx_toa] = fcd.multilateration_toa(deltas_mean, deltas_var, act_bss);
    
    z_mean = deltas_mean;
    z_var = deltas_var;
    
    % ______________________________________________
    z_mean_hist(act_bss,t_idx) = z_mean;
    z_var_hist(act_bss,t_idx) = z_var;
    xy_toa_hist(:,t_idx) = xy_toa;
    xyd_toa_hist(:,t_idx,bss_idx_toa) = xyd_toa;
end
% -------------------------------------------------------------------------
% 4.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist,[0 90 0 90]);

fig2 = figure('Position',[2494 551 560 420]);
show_bss_xy_dist_toa(fig2,target,xyd_toa_hist,[0 90 0 90]);

fig3 = figure('Position',[1932 40 560 420]);
show_measurements(fig3,target.t_vect, z_mean_hist, z_var_hist);

