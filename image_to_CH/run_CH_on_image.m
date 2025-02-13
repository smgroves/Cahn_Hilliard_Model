indir= cd;
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/output";

m = 8;
GridSize = 256;
ny = GridSize;
h = 1/GridSize;
% epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
dt = 1e-5;
max_it = 2000;
boundary = 'neumann';
init_file = sprintf("%s/input/MCF10A_T6I8_chromosome_phi_IC_256.csv",indir);
% init_file = sprintf("%s/input/figure_5_chromosome_phi_IC.csv",indir);

phi0 = readmatrix(init_file);
print_phi = true;
dt_out = 10;
domain_size = 4.4/3.2;
% domain_size = 6.6/3.2;
% grid_ratio = GridSize/128;
epsilon = 0.0067;
% *grid_ratio/domain_size;

domain = [domain_size 0 domain_size 0];

cd ../CHsolvers_package_copy_0206/CahnHilliard_MATLAB_solvers;
pathname = sprintf("%s/MCF10A_SAV_MATLAB_%d_dt_%.2e_Nx_%d_eps_%.5f_",outdir,max_it,dt, GridSize, epsilon);
fprintf("Running SAV solver with parameters: %s\n", pathname);
tStart_SAV = tic;
[t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV_SMG(phi0,...
                                    t_iter = max_it,...
                                    dt = dt,...
                                    m = m,...
                                    boundary = boundary,...
                                    printphi=print_phi,...
                                    pathname=pathname,...
                                    dt_out = dt_out,...
                                    domain = domain);
elapsedTime = toc(tStart_SAV);

fid = fopen('../Job_specs.txt', 'a+');
v = [string(datetime) "FD_spinodal_decomp_smoothed_print" "MATLAB" "FD" GridSize epsilon dt 'NaN' max_it 'NaN' elapsedTime, pathname];
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', v);
fclose(fid);

% writematrix(phi_t(:,:,end),sprintf('%sfinal_phi.csv', pathname));
% writematrix(delta_mass_t,sprintf('%smass.csv', pathname));
% writematrix(E_t,sprintf('%senergy.csv', pathname));

fprintf("Creating movie\n");
filename = strcat(pathname, "movie");
if print_phi
    ch_movie_from_file(strcat(pathname,"phi.csv"), t_out, ny,filename = filename)
else
    ch_movie(phi_t,t_out, filename = filename);
end

% cd ../CHsolvers_package_copy_0206/CahnHilliard_MATLAB_solvers;
% ch_movie_from_file("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/output/HeLa_SAV_MATLAB_2000_dt_1.00e-06_Nx_256_eps_0.00670_phi.csv",...
%                     (0:1:2000)*1e-6*10,...
%                     256,...
%                     filename = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/output/HeLa_SAV_MATLAB_2000_dt_1.00e-06_Nx_256_eps_0.00670_movie")

% cd(indir);



% ch_movie_from_file("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/output/MCF10A_T6I8_SAV_MATLAB_5000_dt_1.00e-06_Nx_128_eps_0.01382_phi.csv",...
%                     (0:1:2000)*1e-5*10,...
%                     128,...
%                     filename = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/output/MCF10A_T6I8_SAV_MATLAB_5000_dt_1.00e-06_Nx_128_eps_0.01382_movie")
