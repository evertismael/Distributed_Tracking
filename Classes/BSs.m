classdef BSs
   properties
      bx              % positions of BSs
      sigma           % ToA noise of furthest
      SNR_db          % Bsbnd noise of furthest to center position
      c = 2.998e8;    % speed of light
      
      pilot_tx
      pilot_rx
      vars
      active_bs       % BS that are active: detected r 
      pilot_grid_r
      
   end
   methods
      function obj = BSs(sigma,SNR_db)
          scene = Params.get_scene();
          obj.bx = scene.bx;
          obj.sigma = sigma;
          obj.SNR_db = SNR_db;
      end
      
      % v1: just distances and gaussian noise
      function [toas, deltas_mean, deltas_var] = compute_toa(obj, x_target, noise_type)
          x_tx = [x_target(1),x_target(3)].';
          d = sqrt(sum((x_tx - obj.bx).^2,1));
          
          t0 = 0;
          if strcmp(noise_type,'same')
              noise = obj.sigma * randn(size(d));
              deltas_var = ((obj.sigma*obj.c)^2)*ones(size(d));
              
          elseif strcmp(noise_type,'sigma_min')
              % compute sigmas based on distance:
              noise = obj.sigma * randn(size(d));
          else
              noise = zeros(size(d));
          end
          
          % final reading:
          toas = d/obj.c + t0 + noise;
          deltas_mean = toas*obj.c;
          
          toas = toas.';
          deltas_mean = deltas_mean.';
          deltas_var = deltas_var.';
          '';
      end
      
      % v2: computes toa from baseband signals
      function obj = gen_pilot_tx(obj)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          % create Tx signal:
          bitstream = randi(2, comm.Nbps * comm.N_pilot, 1) - 1;
          obj.pilot_tx = mapping(bitstream,comm.Nbps,'qam');
          obj.pilot_tx = reshape(obj.pilot_tx,1,1,1,1,comm.N_pilot);
          
          % create signal grid:
          obj.pilot_grid_r = obj.pilot_tx.*exp(1j*comm.phi_rng.*(gs.r));
      end
      function obj = channel_propagation(obj,x_target, noise_type)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          
          x_tx = [x_target(1),x_target(3)].';
          d = sqrt(sum((x_tx - obj.bx).^2,1));
          delta0 = 0; % range offset.
          delta = d + delta0;
          
          % Propagation Channel - Rx signal:
          delta = reshape(delta,1,1,1,length(delta));
          obj.pilot_rx = obj.pilot_tx.*exp(1j*comm.phi_rng.*(delta));
          
          var_n = 10^(-obj.SNR_db/10);
          if strcmp(noise_type,'same')
              w = sqrt(var_n/2)*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));    
              obj.vars = var_n;
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          elseif strcmp(noise_type,'SNR_center')
              % obj.SNR_db: SNR from center to furthest BS
              % distance to center:
              xy_c = [(gs.x_max-gs.x_min)/2, (gs.y_max-gs.y_min)/2].';
              d_c = sqrt(sum((xy_c - obj.bx).^2,1));
              d_c = max(d_c,[],'all');
              % compute SNRs of BSs:
              SNR_lin = 10^(obj.SNR_db/10);
              SNR_i = SNR_lin.*((d_c./d).^2);
              obj.vars = reshape(1./SNR_i,1,1,1,size(d,2));
              w = sqrt(obj.vars/2).*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          elseif strcmp(noise_type,'SNR_20m')
              % obj.SNR_db: SNR when Tx is 20m away from the BS
              d_c = 20;
              % compute SNRs of BSs:
              SNR_lin = 10^(obj.SNR_db/10);
              SNR_i = SNR_lin.*((d_c./d).^2);
              obj.vars = reshape(1./SNR_i,1,1,1,size(d,2));
              w = sqrt(obj.vars/2).*(randn(size(obj.pilot_rx)) + 1j*randn(size(obj.pilot_rx)));
              
              SNR_i_db = 10*log10(SNR_i);
              obj.active_bs = SNR_i_db > -15; % it detects if SNR bigger than 15 db
          else
              w = zeros(size(obj.pilot_rx));
              obj.active_bs = ones(1,size(obj.bx,2)); % TODO
          end
          obj.pilot_rx = obj.pilot_rx + w;
          
      end
      function [toas, deltas_mean, deltas_var, act_bss] = compute_bsbnd_toa(obj)
          
          gs = Params.get_grid_search();
          
          act_bss = find(obj.active_bs == 1); % active bss
          % compute ToA:
          sig_diff = obj.pilot_rx(:,:,:,act_bss,:) - obj.pilot_grid_r;
          llk_delta_bs = (-0.5./obj.vars(:,:,:,act_bss)).*(sum(conj(sig_diff).*sig_diff,5));
          
          % toa_pdf:
          lk = exp(llk_delta_bs - max(llk_delta_bs,[],[1,2,3]));
          llk_norm = trapz(lk,1)*gs.dr;
          
          deltas_mean = (1/llk_norm).*trapz(gs.r.*lk,1)*gs.dr;
          deltas_var = (1/llk_norm).*trapz(((gs.r-deltas_mean).^2).*lk,1)*gs.dr;
          % check min var:
          deltas_var(deltas_var<gs.dr^2) = gs.dr^2;
          
          % output:
          deltas_mean = squeeze(deltas_mean);
          deltas_var = squeeze(deltas_var);
          toas = deltas_mean./obj.c;
          '';
      end
      function [toas, deltas_mean, deltas_var] = refine_toa(obj, prior_mean, prior_var, act_bss)
          
          gs = Params.get_grid_search();
          
          prior = (-0.5./prior_var).*((gs.r - prior_mean).^2);
          prior = reshape(prior, size(prior,1),1,1,size(prior,2));
          % compute ToA:
          sig_diff = obj.pilot_rx(:,:,:,act_bss,:) - obj.pilot_grid_r;
          llk_delta_bs = (-0.5./obj.vars(:,:,:,act_bss)).*(sum(conj(sig_diff).*sig_diff,5));
          llk_delta_bs = llk_delta_bs + prior;
          
          
          % toa_pdf:
          lk = exp(llk_delta_bs - max(llk_delta_bs,[],[1,2,3]));
          llk_norm = trapz(lk,1)*gs.dr;
          
          deltas_mean = (1/llk_norm).*trapz(gs.r.*lk,1)*gs.dr;
          deltas_var = (1/llk_norm).*trapz(((gs.r-deltas_mean).^2).*lk,1)*gs.dr;
          % check min var:
          deltas_var(deltas_var<gs.dr^2) = gs.dr^2;
          
          % output:
          deltas_mean = squeeze(deltas_mean);
          deltas_var = squeeze(deltas_var);
          toas = deltas_mean./obj.c;
      end
      function [toas, deltas_mean, deltas_var] = refine_dist_toa(obj,prior_mean, prior_var, bs_to, act_bss)
          gs = Params.get_grid_search();
          
          % compute ToA:
          sig_diff = obj.pilot_rx(:,:,:,act_bss,:) - obj.pilot_grid_r;
          llk_delta_bs = (-0.5./obj.vars(:,:,:,act_bss)).*(sum(conj(sig_diff).*sig_diff,5));
          
          
          % sum priors:
          prior = zeros(size(llk_delta_bs));
          for tmp_idx = 1: size(act_bss,2)
              bs_idx = act_bss(tmp_idx);
              % collect all priors:
              bs_to_idx = bs_to==bs_idx;
              m_prior_mean = prior_mean(bs_to_idx);
              m_prior_var = prior_var(bs_to_idx);
              
              prior_i = (-0.5./m_prior_var).*((gs.r - m_prior_mean).^2);
              prior_i = reshape(prior_i, size(prior_i,1),1,1,size(prior_i,2));
              prior(:,:,:,tmp_idx) = sum(prior_i,4);
              '';
          end
          
          % add priors:
          llk_delta_bs = llk_delta_bs + prior;
          
          
          % toa_pdf:
          lk = exp(llk_delta_bs - max(llk_delta_bs,[],[1,2,3]));
          llk_norm = trapz(lk,1)*gs.dr;
          
          deltas_mean = (1/llk_norm).*trapz(gs.r.*lk,1)*gs.dr;
          deltas_var = (1/llk_norm).*trapz(((gs.r-deltas_mean).^2).*lk,1)*gs.dr;
          % check min var:
          deltas_var(deltas_var<gs.dr^2) = gs.dr^2;
          
          % output:
          deltas_mean = squeeze(deltas_mean);
          deltas_var = squeeze(deltas_var);
          toas = deltas_mean./obj.c;
          
          
      end
   end
end