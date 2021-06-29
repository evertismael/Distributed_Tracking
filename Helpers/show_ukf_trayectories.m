function show_ukf_trayectories(fig,x_ukf_hist)
figure(fig);
scene = Params.get_scene();
% plot xy-planes for all bss:
for bs_idx = 1:scene.N_bs
    subplot(4,4,bs_idx); hold on;
    ukf_x = x_ukf_hist(1,:,bs_idx);
    ukf_y = x_ukf_hist(3,:,bs_idx);
    t_idx = ukf_x~=0 & ukf_y ~=0;
    
    plot(ukf_x(:,t_idx,:),ukf_y(:,t_idx,:),'r','LineWidth',2);
end

end