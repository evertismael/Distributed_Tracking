clear; clc; close all;
addpath('Classes')  
addpath('Helpers') 
addpath('Targets') 
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
N_iter = 5;


% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -11;


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
    [~, deltas_mean, deltas_var,act_bss] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
    
    % IPE iter 1 onwards:
    xy_iter = xy_toa;
    varxy_iter = varxy_toa;
    for iter_idx = 1:N_iter
        % compute prior:
        [prior_mean, prior_var] = fc.prior_toa(xy_iter, varxy_iter,act_bss);
        
        % refine xy_toa:
        [~, deltas_mean, deltas_var] = bss.refine_toa(prior_mean, prior_var,act_bss);
        
        % re-do multilateration:
        [xy_iter, varxy_iter] = fc.multilateration_toa(deltas_mean, deltas_var,act_bss);
        
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
fig4 = figure('Position',[1943 40 570 413]);
fig5 = figure('Position',[2573 50 570 413]);
fig6 = figure('Position',[3193 52 570 413]);

figs = [fig2,fig3,fig4,fig5,fig6];
show_target_toa_meas(fig1,target,xy_toa_hist,[0 50 0 50]);
for iter_idx = 1:N_iter
    show_target_toa_meas(figs(iter_idx),target,xy_iter_hist(:,:,iter_idx),[0 50 0 50]);
end
''
