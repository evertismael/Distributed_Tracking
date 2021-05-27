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
dt = (1e-2);
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);

% sigma is for toa/roa
% var_n for baseband signal.
SNR_db = -9;


% define BSs and FC
bss = BSs(Inf,SNR_db);
bss = bss.gen_pilot_tx();
fc = FC();

deltas_hist = zeros(scene.N_bs,N_t);
xy_toa_hist = zeros(2,N_t);
for t_idx = 1:N_t
    target.history(:,t_idx);
    
    % noise_type: SNR_center / same
    bss = bss.channel_propagation(target.history(:,t_idx), 'SNR_center');
    
    [~, deltas_mean, deltas_var] = bss.compute_bsbnd_toa();
    
    % compute xy based on toa's at each time-iteration
    [xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var);
    
       
    
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    xy_toa_hist(:,t_idx) = xy_toa;
end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig = figure('Position',[1921 314 1920 988]);
show_target_toa_meas(fig,target,xy_toa_hist);
''


% -------------------------------------------------------------------------
% Extra local Functions
% -------------------------------------------------------------------------

function target = define_target()
x0 = [10,1,5,2].';
target = Mobile(x0);

% add non overlapping (in time) trayectories.
t1 = 0; t2 = 20;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + vx0*(t-t0);
law_y = @(y0,vy0,t0,t) y0 + vy0*(t-t0);
law_vx = @(x0,vx0,t0,t) vx0 + 0*t;
law_vy = @(y0,vy0,t0,t) vy0 + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = 20; t2 = 30;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + vx0*(t-t0);
law_y = @(y0,vy0,t0,t) y0 + 0*t;
law_vx = @(x0,vx0,t0,t) vx0 + 0*t;
law_vy = @(y0,vy0,t0,t) 0 + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = 30; t2 = 50;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + 0*t;
law_y = @(y0,vy0,t0,t) y0 - x0(4)*(t-t0);
law_vx = @(x0,vx0,t0,t) 0 + 0*t;
law_vy = @(y0,vy0,t0,t) -x0(4) + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);
end