classdef mUKF_dist
   properties
      x_est
      phi_est
      P_est
      x_pred
      P_pred
      eig_P_est
      eig_P_pred
      
      f
      A
      Q
      
      h
      R
      
      alpha
      beta
      kappa
      sig_wc
      sig_wm
      
      active
   end
   methods
      function obj = mUKF_dist(alpha, beta, kappa)
          % define f as constant velocity model:
          obj.A = @(T) [1 T 0 0;
                        0 1 0 0;
                        0 0 1 T;
                        0 0 0 1];
          obj.f = @(T,x)obj.A(T)*x;
          Qx = @(T,var_v) [0.5*T^2; T]*var_v*([0.5*T^2; T].');
          obj.Q = @(T,var_v)[Qx(T,var_v) zeros(2,2); zeros(2,2) Qx(T,var_v)];
          
          % define h for range measurements:
          obj.R = @(deltas_var) diag(deltas_var);
          A1 = [1 0 0 0; 0 0 1 0];
          obj.h = @(bx, x_vect) sqrt(sum((A1*x_vect - bx).^2,1)); 
          
          % get sigma weights
          obj.alpha = alpha;
          obj.beta = beta;
          obj.kappa = kappa;
          n = 4;
          [obj.sig_wm, obj.sig_wc] = sigma_weights(n,alpha, beta, kappa);
          '';
          
          obj.active = false;
      end
      
      function obj = correct(obj, bx,z_mean,z_var)
          % sigma points
          n = size(obj.x_est,1);
          sig_x = sigma_points(obj.alpha, obj.kappa, obj.x_pred, obj.P_pred);
                    
          % compute z_ut mean:
          sig_z = zeros(size(bx,2),2*n+1);
          for sig_idx = 1:2*n+1
              sig_z(:,sig_idx) = obj.h(bx, sig_x(:,sig_idx));
          end
          z_ut = sum(obj.sig_wm.*sig_z,2);
          
          % compute z_ut var:
          Pz_ut = zeros(size(sig_z,1),size(sig_z,1));
          for sig_idx = 1:2*n+1
             dz = (sig_z(:,sig_idx) - z_ut);
             Pz_ut = Pz_ut + obj.sig_wc(:,sig_idx)*(dz*dz.');
          end
          Pz_ut = Pz_ut + obj.R(z_var);
          '';
          % compute xz_ut_var
          Pxz_ut = zeros(size(sig_x,1),size(sig_z,1));
          for sig_idx = 1:2*n+1
             dx = (sig_x(:,sig_idx) - obj.x_pred);
             dz = (sig_z(:,sig_idx) - z_ut);
             Pxz_ut = Pxz_ut + obj.sig_wc(:,sig_idx)*(dx*dz.');
          end
          
          
          % kalman gain:
          K = Pxz_ut*pinv(Pz_ut);
          
          % correction:
          obj.phi_est = obj.x_pred + K*(z_mean - z_ut);
          obj.P_est = obj.P_pred - K*Pz_ut*(K.');
          obj.eig_P_est = eig(obj.P_est);
      end
      
      
      function obj = set_phi0(obj, bx, z_means, z_vars)
          % do toa localization with the provided info:
          N_bs = size(bx,2);
          gs = Params.get_grid_search();
          tmp_x = gs.x-reshape(bx(1,:),1,1,N_bs);
          tmp_y = gs.y-reshape(bx(2,:),1,1,N_bs);
          
          gsr = sqrt(tmp_x.^2 + tmp_y.^2);
          
          % assume a time synchronized system => t0 = 0;
          % compute position based on gird-search
          z_means = reshape(z_means,1,1,N_bs);
          z_vars = reshape(z_vars,1,1,N_bs);
          
          llk = (-0.5./z_vars).*(gsr-z_means).^2;
          sum_llk = sum(llk,3);
          
          % get pdf:
          lk = exp(sum_llk-max(sum_llk,[],'all'));
          lk_y = trapz(lk,1)*gs.dx;
          lk_x = trapz(lk,2)*gs.dy;
          pdf_norm = trapz(lk_y,2)*gs.dy;
          
          % compute mean and var
          mean_x = (1/pdf_norm)*trapz(gs.x.*lk_x,1)*gs.dx;
          mean_y = (1/pdf_norm)*trapz(gs.y.*lk_y,2)*gs.dy;
          xy_toa = [mean_x,mean_y].';
          
          tmp_y = trapz(((gs.x - mean_x)).*lk,1)*gs.dx; % for easy computation
          var_x = (1/pdf_norm)*trapz((gs.x - mean_x).^2.*lk_x,1)*gs.dx;
          var_y = (1/pdf_norm)*trapz((gs.y - mean_y).^2.*lk_y,2)*gs.dy;
          var_xy = (1/pdf_norm)*trapz((gs.y - mean_y).*tmp_y,2)*gs.dy;
          varxy_toa = [var_x var_xy; var_xy var_y];  
          
          % check for small variances:
          varxy_toa(varxy_toa < gs.dx^2) = gs.dx^2;
          
          obj.phi_est = ([1 0 0 0; 0 0 1 0].')*xy_toa;
          obj.phi_est = obj.phi_est + [0,rand(1),0,rand(1)].';
          
          % TODO: use varxy_toa to initiate P
          obj.P_est = 10*eye(4);
          
          obj.eig_P_est = eig(obj.P_est);
          obj.eig_P_pred = obj.eig_P_est;
          '';
          
          obj.active = true;
      end
      function obj = predict(obj, T)
          obj.x_pred = obj.f(T,obj.x_est);
          var_v = 0.4;
          obj.P_pred = obj.A(T)*obj.P_est*(obj.A(T).') + obj.Q(T,var_v);
          obj.eig_P_pred = eig(obj.P_pred);
      end
      function obj = mixing_diffusion(obj, mix_type, phis, cs)
            if strcmp(mix_type,'no_mix')
                obj.x_est = obj.phi_est;
            else
                % mix with the Bs in the neighborhood:
                % normalize cs:
                weights = cs./(sum(cs,2));
                obj.x_est = sum(weights.*phis,2);
                '';
            end
      end
   end
end