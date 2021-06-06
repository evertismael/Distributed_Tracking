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
dt = 0.2;%1e-1;
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig = figure('Position',[1921 566 560 420]);
show_target(fig,target,[0 90 0 90]);

fig1 = figure('Position', [3256 346 560 420]);
animate_target(fig1,target,[0 90 0 90]);
