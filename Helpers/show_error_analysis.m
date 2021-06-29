function [mse_i, serr_i] = show_error_analysis(fig,target,x_ukf_hist,active_hist)
figure(fig);
scene = Params.get_scene();

% get distance Error:
subplot(2,1,1);
pos_error = sqrt(sum((target.history - x_ukf_hist).^2,1));
pos_error = squeeze(pos_error);
pos_error(~active_hist.') = 0;
plot(target.t_vect, pos_error);
ylim([0 100]); grid on;
title('\sqrt(distance Error)');

subplot(2,1,2);
similarity_error = zeros(size(pos_error));
for bs_idx = 1:scene.N_bs
    ref_bs = x_ukf_hist(:,:,bs_idx);
    sim_error_i = x_ukf_hist(:,:,1:end~=bs_idx) - ref_bs;
    
    tmp_act = reshape(active_hist(1:end~=bs_idx,:).',1,size(active_hist,2),size(active_hist,1)-1);
    tmp_idx = repmat(tmp_act,4,1,1);
    sim_error_i(~tmp_idx) = 0; 
    
    my_tmp_idx = reshape(active_hist(bs_idx,:).',1,size(active_hist(bs_idx,:),2),size(active_hist(bs_idx,:),1));
    my_tmp_idx = repmat(my_tmp_idx,4,1,11);
    sim_error_i(~my_tmp_idx) = 0; 
    
    
    similarity_error(:,bs_idx) = sqrt(sum(sim_error_i.^2,[1,3]));
end

plot(target.t_vect, similarity_error);
ylim([0 100]); grid on;
title('\sqrt(distance Error)');

mse_i = pos_error.^2;
serr_i = similarity_error.^2;

end