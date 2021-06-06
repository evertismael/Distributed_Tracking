
function target = car_up_right()
x0 = [48.5,0,0,25].';
target = Mobile(x0);
syms t

% TRAYECTORY 1: up 
t1 = 0; t2 = 2.25;
a = 8; b = 0.7; c = .8;
law_vx = @(x_init,v_init,t_init,t) v_init(1) +  0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) - b./(1+exp(-a.*(t-c)))*v_init(2);
law_x = @(x_init,v_init,t_init,t) x_init(1) + 0*t;
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 2: turning to right to start the roundabout.
w = 1.1;
t1 = t2; t2 = t2 + 1.5;
law_vx = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init) - pi/2)+1);
law_vy = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init))+1);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
% 
% TRAYECTORY 3: turning to left following the roundabout.
t1 = t2; t2 = t2 + 5.1;
w = 0.65;
phi = pi/4;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init)+phi));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init)+phi));
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);

% TRAYECTORY 4: turning to right to end the roundabout.
t1 = t2; t2 = t2 + 1;
w = 3;
law_vx = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(cos(w*(t-t_init) + pi/20)-1.8);
law_vy = @(x_init,v_init,t_init,t) 0.5*sqrt(sum(v_init.^2,1))*(sin(w*(t-t_init) - pi/2)-1);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
% 
%TRAYECTORY 5: getting away
t1 = t2; t2 = t2 + 2.5;
a = 3; b = 2; c = t1+1.5;
law_vx = @(x_init,v_init,t_init,t) v_init(1) + b./(1+exp(-a.*(t-c)))*v_init(1);
law_vy = @(x_init,v_init,t_init,t) v_init(2) + 0*t;
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
end