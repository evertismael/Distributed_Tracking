function show_target_and_tracker(fig, target, track, eig_P_est_hist, eig_P_pred_hist)
figure(fig)
    subplot(2,3,1); hold on;
    plot(target.t_vect,target.history(1,:),'DisplayName','true')
    plot(target.t_vect,track(1,:),'DisplayName','ukf')
    title('x'); legend();

    subplot(2,3,4); hold on;
    plot(target.t_vect,target.history(2,:),'DisplayName','true')
    plot(target.t_vect,track(2,:),'DisplayName','ukf')
    title('vx'); ylim([-5,5]); legend();

    subplot(2,3,2); hold on;
    plot(target.t_vect,target.history(3,:),'DisplayName','true')
    plot(target.t_vect,track(3,:),'DisplayName','ukf')
    title('y');legend();

    subplot(2,3,5);hold on;
    plot(target.t_vect,target.history(4,:),'DisplayName','true')
    plot(target.t_vect,track(4,:),'DisplayName','ukf');
    title('vy'); ylim([-5,5]); legend();

    subplot(2,3,6); hold on;
    plot(target.history(1,:),target.history(3,:),'DisplayName','true')
    plot(track(1,:),track(3,:),'DisplayName','ukf')
    
    title('xy-plane');
    xlim([-5,55]); ylim([-5,55]);legend();
    
    
    
    
    
    subplot(2,3,3); hold on;
    plot(target.t_vect,eig_P_est_hist,'--','DisplayName','est')
    plot(target.t_vect,eig_P_pred_hist,'DisplayName','pred')
    title('eig(P)'); legend();
    
end