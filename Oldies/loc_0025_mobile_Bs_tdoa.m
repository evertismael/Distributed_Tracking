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
dt = 1e-1;
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- Define BSs and get ToA measurements from BSs:
% -------------------------------------------------------------------------
scene = Params.get_scene();
N_t = size(target.t_vect,2);

sigma = 20e-9;  % ns; (1/B = 25 ns);
%sigma = 2e-19;

% define BSs and FC
bss = BSs(sigma,0);
fc = FC();

deltas_hist = zeros(scene.N_bs,N_t);
%xy_toa_hist = zeros(2,N_t);
xy_tdoa_hist = zeros(2,N_t);

for t_idx = 1:N_t
    target.history(:,t_idx);
    [~, deltas_mean, deltas_var] = bss.compute_toa(target.history(:,t_idx), 'same');
    
    % compute xy based on toa's at each iteration
    %[xy_toa, varxy_toa] = fc.multilateration_toa(deltas_mean, deltas_var);
    [xy_tdoa, varxy_tdoa] = fc.multilateration_tdoa(deltas_mean, deltas_var);
    
    % ______________________________________________
    % save for history:
    deltas_hist(:,t_idx) = deltas_mean;
    %xy_toa_hist(:,t_idx) = xy_toa;
    xy_tdoa_hist(:,t_idx) = xy_tdoa;
end


% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
%fig1 = figure;
%show_target_toa_meas(fig1,target,xy_toa_hist);
fig = figure;
show_target_toa_meas(fig,target,xy_tdoa_hist);
''
