classdef BSs
   properties
      bx              % positions of BSs
      sigma           % ToA noise of furthest
      SNR_db          % Bsbnd noise of furthest to center position
      c = 2.998e8;    % speed of light
   end
   methods
      function obj = BSs(sigma,SNR_db)
          scene = Params.get_scene();
          obj.bx = scene.bx;
          obj.sigma = sigma;
          obj.SNR_db = SNR_db;
      end
      
      % v1: just distances and gaussian noise
      function [toas, deltas_mean, deltas_var]= compute_toa(obj, x_target, noise_type)
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
      function [toas, deltas_mean, deltas_var] = compute_bpass_toa(obj, x_target, noise_type)
          comm = Params.get_communication();
          gs = Params.get_grid_search();
          
          x_tx = [x_target(1),x_target(3)].';
          d = sqrt(sum((x_tx - obj.bx).^2,1));
          delta0 = 0; % range offset.
          delta = d + delta0;
          
          % create Tx signal:
          bitstream = randi(2, comm.Nbps * comm.N_pilot, 1) - 1;
          pilot_symb_tx = mapping(bitstream,comm.Nbps,'qam');
          pilot_symb_tx = reshape(pilot_symb_tx,1,1,1,1,comm.N_pilot);
          
          % Propagation Channel - Rx signal:
          delta = reshape(delta,1,1,1,length(delta));
          pilot_symb_rx = pilot_symb_tx.*exp(1j*comm.phi_rng.*(delta));
          
          
          var_n = 10^(-obj.SNR_db/10);
          if strcmp(noise_type,'same')
              w = sqrt(var_n/2)*(randn(size(pilot_symb_rx)) + 1j*randn(size(pilot_symb_rx)));    
              vars = var_n;
          elseif strcmp(noise_type,'SNR_center')
              % obj.SNR_db: SNR from center to furthest BS
              % distance to center:
              xy_c = [(gs.x_max-gs.x_min)/2, (gs.y_max-gs.y_min)/2].';
              d_c = sqrt(sum((xy_c - obj.bx).^2,1));
              d_c = max(d_c,[],'all');
              % compute SNRs of BSs:
              SNR_lin = 10^(obj.SNR_db/10);
              SNR_i = SNR_lin.*((d_c./d).^2);
              vars = reshape(1./SNR_i,1,1,1,size(d,2));
              w = sqrt(vars/2).*(randn(size(pilot_symb_rx)) + 1j*randn(size(pilot_symb_rx)));
          else
              w = zeros(size(pilot_symb_rx));
          end
          pilot_symb_rx = pilot_symb_rx + w;
          
          
          % compute ToA:
          pilot_grid_r = pilot_symb_tx.*exp(1j*comm.phi_rng.*(gs.r));
          sig_diff = pilot_symb_rx - pilot_grid_r;
          llk_delta_bs = (-0.5./vars).*(sum(conj(sig_diff).*sig_diff,5));
          
          % toa_pdf:
          lk = exp(llk_delta_bs - max(llk_delta_bs,[],[1,2,3]));
          llk_norm = trapz(lk,1)*gs.dr;
          
          deltas_mean = (1/llk_norm).*trapz(gs.r.*lk,1)*gs.dr;
          deltas_var = (1/llk_norm).*trapz(((gs.r-deltas_mean).^2).*lk,1)*gs.dr;
          
          % output:
          deltas_mean = squeeze(deltas_mean);
          deltas_var = squeeze(deltas_var);
          toas = deltas_mean./obj.c;
      end
   end
end