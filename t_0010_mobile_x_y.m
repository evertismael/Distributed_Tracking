clear; clc; close all;
addpath('Classes')  

% -------------------------------------------------------------------------
% 1.- Define target trajetories
% -------------------------------------------------------------------------
x0 = [0,1,0,2].';
target = Mobile(x0);

% add non overlapping (in time) trayectories.
t1 = 0; t2 = 5;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + vx0*(t-t0);
law_y = @(y0,vy0,t0,t) y0 + vy0*(t-t0);
law_vx = @(x0,vx0,t0,t) vx0 + 0*t;
law_vy = @(y0,vy0,t0,t) vy0 + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = 5; t2 = 10;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + vx0*(t-t0);
law_y = @(y0,vy0,t0,t) y0 + 0*t;
law_vx = @(x0,vx0,t0,t) vx0 + 0*t;
law_vy = @(y0,vy0,t0,t) 0 + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = 10; t2 = 15;
% constant velocity model: but we could have a non-linear system here
law_x = @(x0,vx0,t0,t) x0 + 0*t;
law_y = @(y0,vy0,t0,t) y0 - x0(4)*(t-t0);
law_vx = @(x0,vx0,t0,t) 0 + 0*t;
law_vy = @(y0,vy0,t0,t) -x0(4) + 0*t;
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

% -------------------------------------------------------------------------
% 2.- Generate samples of the trayectory:
% -------------------------------------------------------------------------
dt = 1e-1;
target = target.gen_trayectory(dt);

% -------------------------------------------------------------------------
% 3.- It displays the trayectory - Target only
% -------------------------------------------------------------------------
fig = figure;
show_target(fig,target);
