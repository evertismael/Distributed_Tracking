function show_bss_xy_dist_toa(fig, target, xy_toa_hist, xy_axis)
figure(fig);
scene = Params.get_scene();
% plot xy-planes for all bss:
for bs_idx = 1:scene.N_bs
    subplot(4,4,bs_idx);
    axis(xy_axis); hold on; grid on;
    scatter(xy_toa_hist(1,:,bs_idx),xy_toa_hist(2,:,bs_idx),'+','MarkerEdgeAlpha',0.2)
    % draw scene:
    draw_bss();
    draw_circle(45,45,5); % roundabout:
    draw_road_borders();
    title(num2str(bs_idx));
end

% plot xy over time:
for bs_idx = 1:scene.N_bs
    %tmp_idx = xy_toa_hist(1,:,bs_idx)~=0 & xy_toa_hist(2,:,bs_idx)~=0;
    subplot(4,4,[13,14]); hold on;
    plot(target.t_vect, xy_toa_hist(1,:,bs_idx));
    ylabel('x [m]')
    xlabel('t [s]'); grid on;
    subplot(4,4,[15,16]); hold on;
    plot(target.t_vect, xy_toa_hist(2,:,bs_idx));
    ylabel('y [m]'); grid on;
    xlabel('t [s]')
end
end
