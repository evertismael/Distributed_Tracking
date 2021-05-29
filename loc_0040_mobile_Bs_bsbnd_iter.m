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
N_iter = 5;


% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -10;


% define BSs and FC
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();

deltas_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
xy_iter_hist = zeros(2,N_t,N_iter);

for t_idx = 1:N_t
    target.history(:,t_idx);
    
    % noise_type: SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), 'SNR_center');
    [~, deltas_mean, deltas_var] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_tdoa(deltas_mean, deltas_var);
    
    % IPE iter 1 onwards:
    xy_iter = xy_toa;
    varxy_iter = varxy_toa;
    for iter_idx = 1:N_iter
        % compute prior:
        [prior_mean, prior_var] = fc.prior_toa(xy_iter, varxy_iter);
        
        % refine xy_toa:
        [~, deltas_mean, deltas_var] = bss.refine_toa(prior_mean, prior_var);
        
        % re-do multilateration:
        [xy_iter, varxy_iter] = fc.multilateration_tdoa(deltas_mean, deltas_var);
        
        % save for history:
        xy_iter_hist(:,t_idx,iter_idx) = xy_iter;
    end
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    xy_toa_hist(:,t_idx) = xy_toa;
end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------

fig1 = figure('Position',[1925 847 560 420]);
fig2 = figure('Position',[2573 843 560 420]);
fig3 = figure('Position',[3198 841 560 420]);
fig4 = figure('Position',[1943 342 570 413]);
fig5 = figure('Position',[2573 350 570 413]);
fig6 = figure('Position',[3193 352 570 413]);

figs = [fig2,fig3,fig4,fig5,fig6];
show_target_toa_meas(fig1,target,xy_toa_hist);
for iter_idx = 1:N_iter
    show_target_toa_meas(figs(iter_idx),target,xy_iter_hist(:,:,iter_idx));
end
''


% -------------------------------------------------------------------------
% Extra local Functions
% -------------------------------------------------------------------------

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