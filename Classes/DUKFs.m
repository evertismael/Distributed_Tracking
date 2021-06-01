classdef DUKFs
   properties
       ukfs
       A
   end
   methods
      function obj = DUKFs(N_bs, alpha, beta, kappa)
          obj.ukfs = repmat(microUKF(0,0,0),N_bs,1);
          for mukf_idx = 1:N_bs
              obj.ukfs(mukf_idx) = microUKF(alpha,beta,kappa);
          end
          scene = Params.get_scene();
          obj.A = scene.A;
      end
      
      function obj = set_x0(obj,bx,z_means,z_vars)
          for mukf_idx = 1:size(obj.ukfs,1)
              Ngbrs = obj.A(mukf_idx,:);
              Ngbrs(mukf_idx)=1; % add itself.
              
              ngbr_idx = find(Ngbrs==1);
              mbx = bx(:,ngbr_idx);
              mz_means = z_means(ngbr_idx);
              mz_vars = z_vars(ngbr_idx);
              obj.ukfs(mukf_idx) = obj.ukfs(mukf_idx).set_x0(mbx,mz_means,mz_vars);
          end
      end
      
      function obj = predict(obj,T)
          for mukf_idx = 1:size(obj.ukfs,1)
              obj.ukfs(mukf_idx) = obj.ukfs(mukf_idx).predict(T);
          end
      end
      
      function obj = correct_diffusion(obj, bx, z_means, z_vars)
          for mukf_idx = 1:size(obj.ukfs,1)
              Ngbrs = obj.A(mukf_idx,:);
              Ngbrs(mukf_idx)=1; % add itself.
              
              ngbr_idx = find(Ngbrs==1);
              mbx = bx(:,ngbr_idx);
              mz_means = z_means(ngbr_idx);
              mz_vars = z_vars(ngbr_idx);
              
              obj.ukfs(mukf_idx) = obj.ukfs(mukf_idx).correct_diffusion(mbx,mz_means,mz_vars);
          end
      end
      
      function obj = mixing_diffusion(obj, mix_type)
          % retrieve phi's and eig_P's:
          N_bs = size(obj.ukfs,1);
          phis = zeros(4,N_bs);
          eig_P_ests = zeros(4,N_bs);
          for mukf_idx = 1:size(obj.ukfs,1)
              phis(:,mukf_idx) = obj.ukfs(mukf_idx).phi_est;
              eig_P_ests(:,mukf_idx) = obj.ukfs(mukf_idx).eig_P_est;
          end
          
          % for each ukf: compute cs and mix phis:
          for mukf_idx = 1:size(obj.ukfs,1)
              Ngbrs = obj.A(mukf_idx,:);
              Ngbrs(mukf_idx)=1; % add itself.
              
              ngbr_idx = find(Ngbrs==1);
              mphis = phis(:,ngbr_idx);
              mcs = 1./eig_P_ests(:,ngbr_idx);
              obj.ukfs(mukf_idx) = obj.ukfs(mukf_idx).mixing_diffusion(mix_type, mphis, mcs);
          end
      end

      
      % ___________________________________________________________________
      % helper:
      function [x_ests,eig_P_ests, eig_P_preds] = data_for_history(obj)
          N_bs = size(obj.ukfs,1);
          x_ests = zeros(4,N_bs);
          eig_P_ests = zeros(4,N_bs);
          eig_P_preds = zeros(4,N_bs);
          for mukf_idx = 1:size(obj.ukfs,1)
              x_ests(:,mukf_idx) = obj.ukfs(mukf_idx).x_est;
              eig_P_ests(:,mukf_idx) = obj.ukfs(mukf_idx).eig_P_est;
              eig_P_preds(:,mukf_idx) = obj.ukfs(mukf_idx).eig_P_pred;
          end
      end
   end
end