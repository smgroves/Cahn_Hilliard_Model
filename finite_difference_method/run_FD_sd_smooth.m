indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/";
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function";
m = 16;
GridSize = 128;
dt = 1e-2/(GridSize^2);
n_relax = 4;
h = 1/GridSize;
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
% epsilon = m/(2 * sqrt(2) * atanh(0.9));
gamma = epsilon^2/h^2;
D = GridSize^2;
total_time = (.1/64^2)*500;
% total_time = 0.004;
max_it = round(total_time / dt);
FileName = sprintf("%s/FD_%d_dt_%.2e_Nx_%d_n_relax_%d_gam_%.2e_D_%g",outdir,max_it,dt, GridSize,n_relax, gamma,D);

% max_it = 4000;
% GridSize = 128;
[final_phi,mass_t,E_t] = spinodal_decomp(D,gamma, ...
                "dt",dt,...
                "GridSize",GridSize,...
                "CaptureMode","standard",...
                "NumIterations",max_it,...
                "FrameSpacing",10,...
                'ImgStyle','true',...
                "InputMatrix",sprintf("%s/initial_phi_%d_smooth_n_relax_%d.csv",indir,GridSize,n_relax),...
                "FileName",FileName,...
                "InputType",'phi',...
                "ConstantColorbar", true,...
                "write_residual",false,...
                "write_phi",false);
writematrix(final_phi,sprintf('%s_final_phi.csv', FileName));
writematrix(mass_t,sprintf('%s_mass.csv', FileName));
writematrix(E_t,sprintf('%s_energy.csv', FileName));

fig = figure('visible', 'off');
image(phi_t(:,:,max_it+1),'CDataMapping','scaled'); colorbar; axis square;
set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
clim([-1, 1]);
colormap(redblue(100));
saveas(gca,sprintf('%s_final_phi.png', FileName))
