function animate_target(fig, target, xy_dim)
    figure(fig)
    title('xy-plane')
    axis(xy_dim); grid on;
    hold on;
    % draw scene:
    draw_bss();
    draw_circle(45,45,5); % roundabout:
    draw_road_borders();
    
    for idx = 1:size(target.history(3,:),2)-1
        plot(target.history(1,idx),target.history(3,idx),'.b');
        pause(0.1);
    end  
end


