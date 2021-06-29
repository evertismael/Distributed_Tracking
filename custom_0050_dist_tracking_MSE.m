clear; clc; close all;
addpath('Classes')  
addpath('Helpers')  
addpath('Targets')  
rng(1)
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



% MSE:
N_mse = 300;
se = zeros(N_t,scene.N_bs,N_mse);
serr = zeros(N_t,scene.N_bs,N_mse);

for mse_idx = 1:N_mse

    % sigma is for toa/roa
    % var_n for baseband signal.
    SNR_db = -5;

    % define BSs, FC and UKF_tracker.
    noise_type = 'SNR_20m'; %  SNR_center / same / SNR_20m
    bss = BSs(Inf,SNR_db);
    bss = bss.gen_pilot_tx();
    fc = FC();% just for reference.
    fc_dist = FC_dist();% distributed ToA estimation.

    alpha = 0.2;
    kappa = 3-4; % 3-n 
    beta = 2;
    mix_type = 'no_mix'; % no_mix/eig_P
    ukfs_dist = UKFs_dist(scene.N_bs, alpha, beta, kappa);

    % history variables:
    deltas_mean_hist = zeros(scene.N_bs,N_t);
    deltas_var_hist = zeros(scene.N_bs,N_t);

    deltas_iter_mean_hist = zeros(scene.N_bs,N_t);
    deltas_iter_var_hist = zeros(scene.N_bs,N_t);

    xy_toa_hist = zeros(2,N_t);
    xy_dist_toa_hist = zeros(2,N_t,scene.N_bs);
    xy_dist_iter_hist = zeros(2,N_t,scene.N_bs);
    active_hist = zeros(scene.N_bs,N_t);

    x_est_ukf_hist = zeros(4,N_t,scene.N_bs);
    eig_P_est_hist = zeros(4,N_t,scene.N_bs);
    eig_P_pred_hist = zeros(4,N_t,scene.N_bs);
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
            % prior -> refineToA -> multilateration
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
        % xy_dist_iter
        % tracking:
        ukfs_dist = ukfs_dist.run(z_mean, z_var, act_bss, mix_type, dt);

        % ______________________________________________
        deltas_mean_hist(act_bss,t_idx) = deltas_mean;
        deltas_var_hist(act_bss,t_idx) = deltas_var;

        deltas_iter_mean_hist(act_bss,t_idx) = deltas_iter_mean;
        deltas_iter_var_hist(act_bss,t_idx) = deltas_iter_var;

        xy_toa_hist(:,t_idx) = xy_toa;
        xy_dist_toa_hist(:,t_idx,bss_dist_toa_idx) = xy_dist_toa;
        xy_dist_iter_hist(:,t_idx,bss_dist_toa_idx) = xy_dist_iter;


        [x_ests,eig_P_ests, eig_P_preds, active_vect] = ukfs_dist.data_for_history(act_bss);
        x_est_ukf_hist(:,t_idx,:) = x_ests;
        eig_P_est_hist(:,t_idx,:) = eig_P_ests;
        eig_P_pred_hist(:,t_idx,:) = eig_P_preds;

        active_hist(:,t_idx) = active_vect;
    end
    
    fig = figure;
    [mse_i, serr_i] = show_error_analysis(fig,target, x_est_ukf_hist,active_hist);
    close(fig);
    %______________________________________________________________________
    se(:,:,mse_idx) = mse_i;
    serr(:,:,mse_idx) = serr_i;
    disp(mse_idx);
end
%%
mse = mean(se,3);
mserr = mean(serr,3);
% plot MSE:
fig1 = figure('Position',[1930 634 421 352]);
subplot(2,1,1);
plot(target.t_vect,sqrt(mse));
ylim([0 100]);
xlabel('t [s]'); ylabel('rmse [m]'); grid on;
title('mse');

subplot(2,1,2);
plot(target.t_vect,sqrt(mserr));
ylim([0 100]);
xlabel('t [s]'); ylabel('mean sim err [m]'); grid on;
title('mean_sim_err');
''