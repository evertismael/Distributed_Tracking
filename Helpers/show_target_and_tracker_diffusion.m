function show_target_and_tracker_diffusion(fig, target, track, axis_vect)
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    plot(target.t_vect,squeeze(track(1,:,:)))
    title('x'); legend();

    subplot(2,3,4); hold on;
    plot(target.t_vect,target.history(2,:),'DisplayName','true')
    plot(target.t_vect,squeeze(track(2,:,:)))
    title('vx'); ylim([-5,5]); legend();

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:),'DisplayName','true')
    plot(target.t_vect,squeeze(track(3,:,:)))
    title('y');legend();

    subplot(2,3,5);hold on;
    plot(target.t_vect,target.history(4,:),'DisplayName','true')
    plot(target.t_vect,squeeze(track(4,:,:)))
    title('vy'); ylim([-5,5]); legend();

    subplot(2,3,6); hold on;
    plot(target.history(1,:),target.history(3,:),'DisplayName','true')
    plot(squeeze(track(1,:,:)),squeeze(track(3,:,:)),'DisplayName','ukf')
    
    title('xy-plane');
    axis(axis_vect);
end