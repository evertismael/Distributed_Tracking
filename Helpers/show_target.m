function show_target(fig, target, xy_dim)
figure(fig)
    subplot(2,3,1);
    plot(target.t_vect,target.history(1,:))
    title('x'); grid; xlabel('t[s]')

    subplot(2,3,4);
    plot(target.t_vect,target.history(2,:))
    title('vx'); grid; xlabel('t[s]')

    subplot(2,3,2);
    plot(target.t_vect,target.history(3,:))
    title('y'); grid; xlabel('t[s]')

    subplot(2,3,5);
    plot(target.t_vect,target.history(4,:))
    title('vy'); grid; xlabel('t[s]')

    subplot(2,3,6);
    hold on;
    plot(target.history(1,:),target.history(3,:),'.b');
    plot(target.history(1,1),target.history(3,1),'or');
    
    
    % draw scene:
    draw_bss();
    draw_circle(45,45,5); % roundabout:
    draw_road_borders();
    
    title('xy-plane')
    axis(xy_dim); grid on;
end



