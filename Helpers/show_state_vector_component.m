function show_state_vector_component(fig,target,x_ukf_hist,component_idx)
figure(fig);
scene = Params.get_scene();
% plot xy-planes for all bss:
for bs_idx = 1:scene.N_bs
    subplot(4,4,bs_idx); hold on;
    plot(target.t_vect, target.history(component_idx,:),'DisplayName','True');
    ukf_cmpnt = x_ukf_hist(component_idx,:,bs_idx);
    plot(target.t_vect,ukf_cmpnt(:,:,:),'LineWidth',2,'DisplayName','ukf');
    legend(); grid on;
end
end