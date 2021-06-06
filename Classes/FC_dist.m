classdef FC_dist
   properties
      bx
      A
      sigma
   end
   methods
      function obj = FC_dist()
          scene = Params.get_scene();
          obj.bx = scene.bx;
          obj.A = scene.A;
      end
      % MULTILATERATION (distributed): as many xy_toa as Bss connected to 2
      % or more active Bss
      function [xy_toa, varxy_toa, bss_idx_toa] = multilateration_toa(obj, z_mean, z_var, act_bss)
          gs = Params.get_grid_search();
          
          xy_toa = double.empty(2,0);
          varxy_toa = double.empty(2,2,0);
          bss_idx_toa = double.empty(1,0);
          
          for tmp_idx = 1:size(act_bss,2) 
              act_bs_idx = act_bss(tmp_idx);
              % get neigbours:
              Ngbrs = obj.A(act_bs_idx,:);
              Ngbrs(act_bs_idx)=1; % add itself.
              ngbr_idx = find(Ngbrs==1);
              act_ngbrs = intersect(act_bss,ngbr_idx);
              
              N_ngbhd = size(act_ngbrs,2);  
              if N_ngbhd >= 3 % enough data -> do localization:
                  tmp_x = gs.x-reshape(obj.bx(1,act_ngbrs),1,1,N_ngbhd);
                  tmp_y = gs.y-reshape(obj.bx(2,act_ngbrs),1,1,N_ngbhd);

                  gsr = sqrt(tmp_x.^2 + tmp_y.^2);
                  
                  % select z from act ngbhd:
                  [~, act_ngbrs_idx] = intersect(act_bss,act_ngbrs);
                  
                  % assume a time synchronized system => t0 = 0;
                  % compute position based on gird-search
                  mz_mean = z_mean(act_ngbrs_idx);
                  mz_var = z_var(act_ngbrs_idx);
                                    
                  mz_mean = reshape(mz_mean,1,1,N_ngbhd);
                  mz_var = reshape(mz_var,1,1,N_ngbhd);
                

                  llk = (-0.5./mz_var).*(gsr-mz_mean).^2;
                  sum_llk = sum(llk,3);
                  
                  % get pdf:
                  lk = exp(sum_llk-max(sum_llk,[],'all'));
                  lk_y = trapz(lk,1)*gs.dx;
                  lk_x = trapz(lk,2)*gs.dy;
                  pdf_norm = trapz(lk_y,2)*gs.dy;

                  % compute mean and var
                  mean_x = (1/pdf_norm)*trapz(gs.x.*lk_x,1)*gs.dx;
                  mean_y = (1/pdf_norm)*trapz(gs.y.*lk_y,2)*gs.dy;
                  xy_toa_bs = [mean_x,mean_y].';

                  tmp_y = trapz(((gs.x - mean_x)).*lk,1)*gs.dx; % for easy computation
                  var_x = (1/pdf_norm)*trapz((gs.x - mean_x).^2.*lk_x,1)*gs.dx;
                  var_y = (1/pdf_norm)*trapz((gs.y - mean_y).^2.*lk_y,2)*gs.dy;
                  var_xy = (1/pdf_norm)*trapz((gs.y - mean_y).*tmp_y,2)*gs.dy;
                  varxy_toa_bs = [var_x var_xy; var_xy var_y];  

                  % check for small variances:
                  varxy_toa_bs(varxy_toa_bs < gs.dx^2) = gs.dx^2;
                  
                  xy_toa = cat(2,xy_toa,xy_toa_bs);
                  varxy_toa = cat(3,varxy_toa,varxy_toa_bs);
                  bss_idx_toa = cat(2,bss_idx_toa,act_bs_idx);
                  '';
              end
          end
          '';
      end
     
      function [prior_mean, prior_var, bs_from, bs_to] = prior_toa(obj, xy_dist_iter, varxy_dist_iter,bss_dist_idx, act_bss)
          prior_mean = double.empty(1,0);
          prior_var = double.empty(1,0);
          bs_from = double.empty(1,0);
          bs_to = double.empty(1,0);
          
          
          % compute prior for each received toa (excluding that toa):
          for tmp_bs_idx = 1:size(bss_dist_idx,2)
              bss_dist = bss_dist_idx(tmp_bs_idx);
              xy_mean = xy_dist_iter(:,tmp_bs_idx);
              xy_var = varxy_dist_iter(:,:,tmp_bs_idx);
              
              % get active neighbors:
              Ngbrs = obj.A(bss_dist,:);
              ngbr_idx = find(Ngbrs==1);
              act_ngbrs = intersect(act_bss,ngbr_idx);
              
              % prior_mean:
              prior_mean_i = sqrt(sum((xy_mean - obj.bx(:,act_ngbrs)).^2,1));
              % prior_var:
              xy_diff = obj.bx(:,act_ngbrs) - xy_mean;
              tmp = [sum(xy_diff.*xy_var(:,1),1); sum(xy_diff.*xy_var(:,2),1)];
              prior_var_i = (1./(prior_mean_i.^2)).*sum(tmp.*xy_diff,1);
            
              % add to the output:
              prior_mean = cat(2,prior_mean,prior_mean_i);
              prior_var = cat(2,prior_var,prior_var_i);
              bs_to = cat(2,bs_to,act_ngbrs);
              bs_from = cat(2,bs_from,bss_dist*ones(size(act_ngbrs)));
          end
          '';
      end
   end
end