classdef Mobile
   properties
      x0
      history
      t_vect
      moves_inter
      moves_law
   end
   methods
      function obj = Mobile(x0)
          obj.x0 = x0;
          obj.moves_inter = double.empty(2, 0);
          obj.moves_law = repmat(Law(0,0,0,0),1,0);
      end
      
      function obj = add_trayectory(obj,t1,t2,law)
          obj.moves_inter(:,end+1) = [t1,t2].';
          obj.moves_law(end+1) = law;
      end
      
      function obj = gen_trayectory(obj,dt)
          obj.t_vect = min(obj.moves_inter,[],'all'):dt:max(obj.moves_inter,[],'all')-dt;
          obj.history = zeros(4,size(obj.t_vect,2));
          
          % generate by moves_law
          for law_idx = 1:size(obj.moves_inter,2)
             t0 = obj.moves_inter(1,law_idx);
             t_idxs = t0 <= obj.t_vect & obj.t_vect < obj.moves_inter(2,law_idx);
             if law_idx == 1
                 x_init = [1,0,0,0;0,0,1,0]*obj.x0;
                 v_init = [0,1,0,0;0,0,0,1]*obj.x0;
             else
                 t_idx_init = find(t_idxs>0,1);
                 x_init = [1,0,0,0;0,0,1,0]*obj.history(:,t_idx_init-1);
                 v_init = [0,1,0,0;0,0,0,1]*obj.history(:,t_idx_init-1);
             end
             
             x = obj.moves_law(law_idx).law_x(x_init,v_init,t0,obj.t_vect(t_idxs));
             y = obj.moves_law(law_idx).law_y(x_init,v_init,t0,obj.t_vect(t_idxs));
             vx = obj.moves_law(law_idx).law_vx(x_init,v_init,t0,obj.t_vect(t_idxs));
             vy = obj.moves_law(law_idx).law_vy(x_init,v_init,t0,obj.t_vect(t_idxs));
             obj.history(:,t_idxs)=[x;vx;y;vy];             
          end 
      end
   end
end