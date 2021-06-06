function target = oval_tragectory()
x0 = [20,2,10,0].';
target = Mobile(x0);
syms t
% add non overlapping (in time) trayectories.
t1 = 0; t2 = 5;
% constant velocity model: but we could have a non-linear system here
law_vx = @(x_init,v_init,t_init,t) v_init(1) + 0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) + 0*t;
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);


% constant velocity model: but we could have a non-linear system here
w0 = 0.15;
t1 = 5; t2 = 5 + pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init));
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init));
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

t1 = t2; t2 = t2 + 5;
% constant velocity model: but we could have a non-linear system here
law_vx = @(x_init,v_init,t_init,t) v_init(1) + 0*t;
law_vy = @(x_init,v_init,t_init,t) v_init(2) + 0*t;
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);

target = target.add_trayectory(t1,t2,law);

% constant velocity model: but we could have a non-linear system here
w0 = 0.15;
t1 = t2; t2 = t2 + pi/w0;
law_vx = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*cos(w0*(t-t_init) + pi);
law_vy = @(x_init,v_init,t_init,t) sqrt(sum(v_init.^2))*sin(w0*(t-t_init) + pi);
law_x = @(x_init,v_init,t_init,t) x_init(1) + int(law_vx(x_init,v_init,t_init,t),t,t_init,t);
law_y = @(x_init,v_init,t_init,t) x_init(2) + int(law_vy(x_init,v_init,t_init,t),t,t_init,t);
law = Law(law_x,law_y,law_vx, law_vy);
target = target.add_trayectory(t1,t2,law);
end