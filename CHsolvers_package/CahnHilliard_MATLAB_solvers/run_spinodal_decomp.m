indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/";
n_relax = 4;
m = 8;
GridSize = 128;
h = 1/GridSize;
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
D = GridSize^2;
% total_time = .2;
dt = 5.5e-6;
% max_it = round(total_time / dt);
max_it = 2000;
boundary = 'neumann';
init_file = sprintf("%s/initial_phi_%d_smooth_n_relax_%d.csv",indir,GridSize, n_relax);
phi0 = readmatrix(init_file);
print_phi = true;
% #################################################
% RUN SAV SOLVER 
% #################################################

% outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/SAV/output/spinodal_smooth_relax_function/";
% pathname = sprintf("%s/SAV_MATLAB_%d_dt_%.2e_Nx_%d_n_relax_%d_",outdir,max_it,dt, GridSize, n_relax);
% tStart_SAV = tic;
% [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV(phi0,...
%                                     t_iter = max_it,...
%                                     dt = dt,...
%                                     m = m,...
%                                     boundary = boundary,...
%                                     printphi=print_phi,...
%                                     pathname=pathname,...
%                                     dt_out = 10);
% elapsedTime = toc(tStart_SAV);

% fid = fopen('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/CHsolvers_package/Job_specs.txt', 'a+');
% v = [string(datetime) "SAV_spinodal_decomp_smoothed_print" "MATLAB" "SAV" GridSize epsilon dt 'NaN' max_it 'NaN' elapsedTime];
% fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', v);
% fclose(fid);

% writematrix(phi_t(:,:,end),sprintf('%sfinal_phi.csv', pathname));
% writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
% writematrix(E_t,sprintf('%senergy.csv', pathname));

% filename = strcat(pathname, "movie");
% ch_movie(phi_t,t_out, filename = filename);


% #################################################
% RUN NMG SOLVER 
% #################################################

% outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/matlab_multigrid/output/spinodal_smooth_relax_function/";
% pathname = sprintf("%s/NMG_MATLAB_%d_dt_%.2e_Nx_%d_n_relax_%d_",outdir,max_it,dt, GridSize, n_relax);
% tStart_NMG = tic;
% [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_NMG(phi0,...
%                                     t_iter = max_it,...
%                                     dt = dt,...
%                                     m = m,...
%                                     boundary = boundary,...
%                                     printphi=print_phi,...
%                                     pathname=pathname,...
%                                     dt_out = 10);
% elapsedTime = toc(tStart_NMG);

% fid = fopen('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/CHsolvers_package/Job_specs.txt', 'a+');
% v = [string(datetime) "NMG_spinodal_decomp_smoothed_print" "MATLAB" "NMG" GridSize epsilon dt 1e-5 max_it 1e4 elapsedTime];
% fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', v);
% fclose(fid);

% writematrix(phi_t(:,:,end),sprintf('%sfinal_phi.csv', pathname));
% writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
% writematrix(E_t,sprintf('%senergy.csv', pathname));

% filename = strcat(pathname, "movie");
% ch_movie(phi_t,t_out, filename = filename);



% % #################################################
% % RUN FD SOLVER 
% % #################################################

% outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function/";
% pathname = sprintf("%s/FD_MATLAB_%d_dt_%.2e_Nx_%d_n_relax_%d_",outdir,max_it,dt, GridSize, n_relax);
% tStart_FD = tic;
% [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_FD(phi0,...
%                                     t_iter = max_it,...
%                                     dt = dt,...
%                                     m = m,...
%                                     boundary = boundary,...
%                                     printphi=print_phi,...
%                                     pathname=pathname,...
%                                     dt_output = 10);
% elapsedTime = toc(tStart_FD);

% fid = fopen('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/CHsolvers_package/Job_specs.txt', 'a+');
% v = [string(datetime) "spinodal_decomp_smoothed_print" "MATLAB" "FD" GridSize epsilon dt 'NaN' max_it 'NaN' elapsedTime];
% fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', v);
% fclose(fid);

% writematrix(phi_t(:,:,end),sprintf('%sfinal_phi.csv', pathname));
% writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
% writematrix(E_t,sprintf('%senergy.csv', pathname));

% filename = strcat(pathname, "movie");
% ch_movie(phi_t,t_out, filename = filename);

