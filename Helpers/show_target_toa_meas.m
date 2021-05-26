function show_target_toa_meas(fig, target, xy_toa_hist)
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    scatter(target.t_vect,xy_toa_hist(1,:),'+','MarkerEdgeAlpha',0.2,'DisplayName','multilat')
    title('x'); legend();

    subplot(2,3,4);
    plot(target.t_vect,target.history(2,:))
    title('vx')

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:));
    scatter(target.t_vect,xy_toa_hist(2,:),'+','MarkerEdgeAlpha',0.2,'DisplayName','multilat')
    title('y')

    subplot(2,3,5);
    plot(target.t_vect,target.history(4,:))
    title('vy')

    subplot(2,3,6);
    hold on;
    plot(target.history(1,:),target.history(3,:),'.b');
    scatter(xy_toa_hist(1,:),xy_toa_hist(2,:),'+','MarkerEdgeAlpha',0.3,'DisplayName','multilat')
    plot(target.history(1,1),target.history(3,1),'or'); % initial point
    
    title('xy-plane')
    xlim([-5,55]); ylim([-5,55]);
end