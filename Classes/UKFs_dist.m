classdef UKFs_dist
   properties
       ukfs
       A
   end
   methods
      function obj = UKFs_dist(N_bs, alpha, beta, kappa)
          obj.ukfs = repmat(mUKF_dist(0,0,0),N_bs,1);
          for ukf_idx = 1:N_bs
              obj.ukfs(ukf_idx) = mUKF_dist(alpha,beta,kappa);
          end
          scene = Params.get_scene();
          obj.A = scene.A;
      end
      
      function obj = run(obj, z_mean, z_var, act_bss, mix_type, dt)
          scene = Params.get_scene();
          
          for tmp_idx = 1:size(act_bss,2)
              ukf_idx = act_bss(tmp_idx);
              % get meas for that bss
              Ngbrs = obj.A(ukf_idx,:);
              Ngbrs(ukf_idx)=1; % add itself.
              ngbr_idx = find(Ngbrs==1);
              [act_ngbrs,  act_ngbrs_idxs] = intersect(act_bss,ngbr_idx);
            
              mz_mean = z_mean(act_ngbrs_idxs);
              mz_var = z_var(act_ngbrs_idxs);
              mbx = scene.bx(:,act_ngbrs);
              
              % only if it has enough data to localize:
              if size(act_ngbrs,2)>=3
                  % correct / set phi_0:
                  if obj.ukfs(ukf_idx).active
                      obj.ukfs(ukf_idx) = obj.ukfs(ukf_idx).correct(mbx,mz_mean,mz_var);
                  else 
                     obj.ukfs(ukf_idx) = obj.ukfs(ukf_idx).set_phi0(mbx,mz_mean,mz_var);
                  end
              end
          end
          
          
          for tmp_idx = 1:size(act_bss,2)
              ukf_idx = act_bss(tmp_idx);
              % get meas for that bss
              Ngbrs = obj.A(ukf_idx,:);
              Ngbrs(ukf_idx)=1; % add itself.
              ngbr_idx = find(Ngbrs==1);
              [act_ngbrs,  ~] = intersect(act_bss,ngbr_idx);
              
              % only if it has enough data to localize:
              if size(act_ngbrs,2)>=3
                  % mix diffusion:
                  N_act_ngbrs = size(act_ngbrs,2);
                  % check active and tracking:
                  act_ngbrs_trk = [];
                  for tmp_ukf_idx = 1:N_act_ngbrs
                      mukf_idx = act_ngbrs(tmp_ukf_idx);
                      if obj.ukfs(mukf_idx).active
                         act_ngbrs_trk = cat(2,act_ngbrs_trk,mukf_idx);
                      end
                  end
                  
                  % use neighbors that are active-tracking
                  '';
                  N_act_ngbrs_trk = size(act_ngbrs_trk,2);
                  mphis = zeros(4,N_act_ngbrs_trk);
                  meig_P_ests = zeros(4,N_act_ngbrs_trk);
                  for tmp_ukf_idx = 1:N_act_ngbrs_trk
                      mukf_idx = act_ngbrs_trk(tmp_ukf_idx);
                      mphis(:,tmp_ukf_idx) = obj.ukfs(mukf_idx).phi_est;
                      meig_P_ests(:,tmp_ukf_idx) = obj.ukfs(mukf_idx).eig_P_est;
                      '';
                  end
                  % do the mixing:
                  mcs = 1./meig_P_ests;
                  obj.ukfs(ukf_idx) = obj.ukfs(ukf_idx).mixing_diffusion(mix_type, mphis, mcs);

                  % predict:
                  obj.ukfs(ukf_idx) = obj.ukfs(ukf_idx).predict(dt);
              end
          end
          '' ;
         
      end
      
      % ___________________________________________________________________
      % helper:
      function [x_ests,eig_P_ests, eig_P_preds, active_vect] = data_for_history(obj, act_bss)
          N_bs = size(obj.ukfs,1);
          x_ests = zeros(4,N_bs);
          eig_P_ests = zeros(4,N_bs);
          eig_P_preds = zeros(4,N_bs);
          active_vect = zeros(1,N_bs);
          for mukf_idx = 1:size(obj.ukfs,1)
              if ~isempty(obj.ukfs(mukf_idx).x_est)
                  x_ests(:,mukf_idx) = obj.ukfs(mukf_idx).x_est;
              end
              if ~isempty(obj.ukfs(mukf_idx).eig_P_est)
                  eig_P_ests(:,mukf_idx) = obj.ukfs(mukf_idx).eig_P_est;
              end
              if ~isempty(obj.ukfs(mukf_idx).eig_P_pred)
                  eig_P_preds(:,mukf_idx) = obj.ukfs(mukf_idx).eig_P_pred;
              end
              
              active_vect(:,mukf_idx) = obj.ukfs(mukf_idx).active && ismember(mukf_idx,act_bss);
          end
      end
   end
end