clear; clc; close all;
addpath('Classes')  
addpath('Helpers')  

% -------------------------------------------------------------------------
% 1.- Define target trajetories
% -------------------------------------------------------------------------
target = define_target();

% -------------------------------------------------------------------------
% 2.- Generate samples of the trayectory:
% -------------------------------------------------------------------------
dt = (1e-1);
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_iter = 1;

% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -13;


% define BSs, FC and UKF_tracker.
noise_type = 'SNR_center'; %  SNR_center / same
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();% just for reference.

alpha = 0.2;
kappa = 3-4; % 3-n 
beta = 2;
N_bs = 4;
dukfs = DUKFs(N_bs, alpha, beta, kappa);

% history variables:
deltas_mean_hist = zeros(scene.N_bs,N_t);
deltas_var_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
xy_ukf_hist = zeros(4,N_t,N_bs);
eig_P_est_hist = zeros(4,N_t,N_bs);
eig_P_pred_hist = zeros(4,N_t,N_bs);


% first sample initiates the tracker:
for t_idx = 1:N_t
    % read sensor
    noise_type = 'SNR_center'; %  SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var);
    
    % iterations:
    xy_iter = xy_toa;
    varxy_iter = varxy_toa;
    for iter_idx = 1:N_iter
        % prior -> refineToA -> multilateration:
        [prior_mean, prior_var] = fc.prior_toa(xy_iter, varxy_iter);
        [~, deltas_iter_mean, deltas_iter_var] = bss.refine_toa(prior_mean, prior_var);
        [xy_iter, varxy_iter] = fc.multilateration_toa(deltas_iter_mean, deltas_iter_var);
    end
    
    if N_iter == 0
        z_mean = deltas_mean;
        z_var = deltas_var;
    else
        z_mean = deltas_iter_mean;
        z_var = deltas_iter_var;
    end
    
    
    if t_idx == 1
        % initiate track:
        dukfs = dukfs.set_x0(bss.bx,z_mean,z_var);
    else
        % dukf: correct with z shared by diffusion.
        dukfs = dukfs.correct_diffusion(bss.bx,z_mean,z_var);
        % mix phi's shared by diffusion:
        mix_type = 'no_mix'; % no_mix/eig_P
        dukfs = dukfs.mixing_diffusion(mix_type);
    end
    % predict: normaly.
    dukfs = dukfs.predict(dt);
    '';
    
    % ______________________________________________
    % save for history:
    % measurements, cenralized loc, ukfilters, eig_P
    deltas_mean_hist(:,t_idx) = deltas_mean;
    deltas_var_hist(:,t_idx) = deltas_var;
    xy_toa_hist(:,t_idx) = xy_toa;
    
    [x_ests,eig_P_ests, eig_P_preds] = dukfs.data_for_history();
    xy_ukf_hist(:,t_idx,:) = x_ests;
    eig_P_est_hist(:,t_idx,:) = eig_P_ests;
    eig_P_pred_hist(:,t_idx,:) = eig_P_preds;
end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------

fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist);
fig2 = figure('Position',[1946 340 560 420]);
show_measurements(fig2, target.t_vect, deltas_mean_hist, deltas_var_hist);
fig3 = figure('Position',[2491 854 570 413]);
show_target_and_tracker_diffusion(fig3, target, xy_ukf_hist);
fig4 = figure('Position',[2506 339 560 420]);
show_covariances_diffusion(fig4, target, eig_P_est_hist,eig_P_pred_hist);

% -------------------------------------------------------------------------
% Extra local Functions
% -------------------------------------------------------------------------

function target = define_target_()
x0 = [20,2,10,0].';
target = Mobile(x0);

% constant velocity model: but we could have a non-linear system here
w0 = 0.15;
t1 = 0; t2 = 2*pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init));
law_x = @(x_init,v_init,t_init,t) x_init(1) + (sqrt(sum(v_init.^2))/w0)*sin(w0*(t-t_init));
law_y = @(x_init,v_init,t_init,t) x_init(2) - (sqrt(sum(v_init.^2))/w0)*(cos(w0*(t-t_init)) - 1);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);
end

function target = define_target()
x0 = [20,2,10,0].';
target = Mobile(x0);

% add non overlapping (in time) trayectories.
t1 = 0; t2 = 5;
% constant velocity model: but we could have a non-linear system here
law_x = @(x_init,v_init,t_init,t) x_init(1) + v_init(1)*(t-t_init);
law_y = @(x_init,v_init,t_init,t) x_init(2) + v_init(2)*(t-t_init);
law_vx = @(x_init,v_init,t_init,t) v_init(1) + 0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);


% constant velocity model: but we could have a non-linear system here
w0 = 0.15;
t1 = 5; t2 = 5 + pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init));
law_x = @(x_init,v_init,t_init,t) x_init(1) + (sqrt(sum(v_init.^2))/w0)*sin(w0*(t-t_init));
law_y = @(x_init,v_init,t_init,t) x_init(2) - (sqrt(sum(v_init.^2))/w0)*(cos(w0*(t-t_init)) - 1);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = t2; t2 = t2 + 5;
% constant velocity model: but we could have a non-linear system here
law_x = @(x_init,v_init,t_init,t) x_init(1) + v_init(1)*(t-t_init);
law_y = @(x_init,v_init,t_init,t) x_init(2) + v_init(2)*(t-t_init);
law_vx = @(x_init,v_init,t_init,t) v_init(1) + 0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

% constant velocity model: but we could have a non-linear system here
w0 = 0.15;
t1 = t2; t2 = t2 + pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init) + pi);
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init) + pi);
law_x = @(x_init,v_init,t_init,t) x_init(1) + (sqrt(sum(v_init.^2))/w0)*sin(w0*(t-t_init) + pi);
law_y = @(x_init,v_init,t_init,t) x_init(2) + (sqrt(sum(v_init.^2))/w0)*(cos(w0*(t-t_init))-1);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
end