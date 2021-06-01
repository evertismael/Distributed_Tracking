function show_covariances_diffusion(fig, target, eig_P_est_hist, eig_P_pred_hist)
    figure(fig)
    for ukf_idx = 1:size(eig_P_est_hist,3)
        subplot(2,2,ukf_idx); hold on;
        plot(target.t_vect,eig_P_est_hist(:,:,ukf_idx),'--','DisplayName','est')
        plot(target.t_vect,eig_P_pred_hist(:,:,ukf_idx),'DisplayName','pred')
        title('eig(P)'); legend();
    end
    
    
end