classdef Params
   methods (Static)
       function scene = get_scene()
           persistent bx;
           persistent network_type;
           
           network_type = 'ring'; % fully/ring
           
           scene = {};
           bx = [0,0; 50,0; 50,50; 0,50].';
           scene.bx = bx;
           scene.N_bs = size(bx,2);
           
           % Adj matrix:
           if strcmp(network_type,'fully')
               % Fully connected netxork:
               scene.bs_from = [1 1 1, 2 2 2, 3 3 3, 4 4 4];
               scene.bs_to  =  [2 3 4, 1 3 4, 1 2 4, 1 2 3];
               scene.A = zeros(scene.N_bs,scene.N_bs);
               idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
               scene.A(idx_from_to) = 1;
           elseif strcmp(network_type,'ring')
               % Ring-Gossip network:
               scene.bs_from = [1 1, 2 2, 3 3, 4 4];
               scene.bs_to  =  [2 4, 1 3, 2 4, 1 3];
               scene.A = zeros(scene.N_bs,scene.N_bs);
               idx_from_to = sub2ind(size(scene.A), scene.bs_from,scene.bs_to);
               scene.A(idx_from_to) = 1;
           else
               error('network_type not recognized');
           end
       end
       
       
       function gs = get_grid_search()
           gs = {};
           gs.x_min = 0; gs.x_max = 50; gs.Nx = 50;
           gs.y_min = 0; gs.y_max = 50; gs.Ny = 50;
           gs.r_min = 0; gs.r_max = 150; gs.Nr = 200;
           
           gs.dx = (gs.x_max - gs.x_min)/gs.Nx;
           gs.dy = (gs.y_max - gs.y_min)/gs.Ny;
           gs.dr = (gs.r_max - gs.r_min)/gs.Nr;
           
           gs.x = (gs.x_min:gs.dx:gs.x_max-gs.dx).';
           gs.y = (gs.y_min:gs.dy:gs.y_max-gs.dy);
           gs.r = (gs.r_min:gs.dr:gs.r_max-gs.dr).';
       end
       
       function comm = get_communication()
           comm = {};
           comm.c = 2.998e8;
           comm.fc = 2e9;
           comm.B = 40e6;
           comm.N_subcrr = 1024;
           comm.pilot_spacing = 16;
           comm.N_pilot = comm.N_subcrr/comm.pilot_spacing;
           comm.deltaPhi = 2*pi*comm.B/(comm.N_pilot*comm.c);
           comm.Nbps = 2;
           
           tmp.pilot_idx_rng = (1:comm.N_pilot);
           tmp.phi_rng = comm.deltaPhi*tmp.pilot_idx_rng;
           comm.phi_rng = reshape(tmp.phi_rng,[1,1,1,1,comm.N_pilot]);
           
           
       end
   end
end