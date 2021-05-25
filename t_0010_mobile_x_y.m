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
figure;
subplot(2,3,1);
plot(target.t_vect,target.history(1,:))
title('x')

subplot(2,3,4);
plot(target.t_vect,target.history(2,:))
title('vx')

subplot(2,3,2);
plot(target.t_vect,target.history(3,:))
title('y')

subplot(2,3,5);
plot(target.t_vect,target.history(4,:))
title('vy')

subplot(2,3,6);
title('xy-plane')
xlim([0,12]);
ylim([0,12]);

hold on;
for t_idx = 1:1:size(target.t_vect,2)
    plot(target.history(1,t_idx),target.history(3,t_idx),'.b');
    pause(0.01)
end
