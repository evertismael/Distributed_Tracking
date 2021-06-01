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
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);
N_iter = 1;

% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -13;


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

deltas_iter_mean_hist = zeros(scene.N_bs,N_t);
deltas_iter_var_hist = zeros(scene.N_bs,N_t);

xy_toa_hist = zeros(2,N_t);
xy_iter_hist = zeros(2,N_t);
xy_ukf_hist = zeros(4,N_t);
eig_P_est_hist = zeros(4,N_t);
eig_P_pred_hist = zeros(4,N_t);


% first sample initiates the tracker:
for t_idx = 1:N_t
    %target.history(:,t_idx)
    noise_type = 'SNR_center'; %  SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), noise_type);
    [~, deltas_mean, deltas_var] = bss.compute_bsbnd_toa();
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var);
    
    % IPE iter 2 onwards:
    xy_iter = xy_toa;
    varxy_iter = varxy_toa;
    for iter_idx = 1:N_iter
        % prior -> refineToA -> multilateration -> saveHist:
        [prior_mean, prior_var] = fc.prior_toa(xy_iter, varxy_iter);
        [~, deltas_iter_mean, deltas_iter_var] = bss.refine_toa(prior_mean, prior_var);
        [xy_iter, varxy_iter] = fc.multilateration_toa(deltas_iter_mean, deltas_iter_var);
        xy_iter_hist(:,t_idx) = xy_iter;
    end
    
    
    if t_idx == 1
        % initialization:
        ukf_fcnp = ukf_fcnp.set_x0(xy_iter);
    else
        if N_iter >= 1
            z_mean = deltas_iter_mean;
            z_var = deltas_iter_var;
        else
            z_mean = deltas_mean;
            z_var = deltas_var;
        end
        % tracking:
        dt = target.t_vect(t_idx) - target.t_vect(t_idx-1); % in case is not continuous sampling
        ukf_fcnp = ukf_fcnp.predict(dt);
        ukf_fcnp = ukf_fcnp.correct(z_mean,z_var);
    end
    
    % ______________________________________________
    % save for history:
    deltas_mean_hist(:,t_idx) = deltas_mean;
    deltas_var_hist(:,t_idx) = deltas_var;
    if N_iter>=1
        deltas_iter_mean_hist(:,t_idx) = deltas_iter_mean;
        deltas_iter_var_hist(:,t_idx) = deltas_iter_var;
    end
    xy_toa_hist(:,t_idx) = xy_toa;
    xy_ukf_hist(:,t_idx) = ukf_fcnp.x_est;
    eig_P_est_hist(:,t_idx) = ukf_fcnp.eig_P_est;
    eig_P_pred_hist(:,t_idx) = ukf_fcnp.eig_P_pred;

end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------

fig1 = figure('Position',[1925 847 560 420]);
show_target_toa_meas(fig1,target,xy_toa_hist);

if N_iter>=1
    fig2 = figure('Position',[2577 858 560 420]);
    show_target_toa_meas(fig2,target,xy_iter_hist);
    fig4 = figure('Position',[2577 346 560 420]);
    show_measurements(fig4, target.t_vect, deltas_iter_mean_hist, deltas_iter_var_hist);
end

fig3 = figure('Position',[1946 340 560 420]);
show_measurements(fig3, target.t_vect, deltas_mean_hist, deltas_var_hist);


fig5 = figure('Position',[3203 855 560 420]);
show_target_and_tracker(fig5, target, xy_ukf_hist,eig_P_est_hist,eig_P_pred_hist);
''


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