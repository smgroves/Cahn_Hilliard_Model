
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/IC/";
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/matlab_multigrid/output/large_and_small_droplets";
r1 = 20;
r2 = 30;
space = 10;
% epsilon = m/(2 * sqrt(2) * atanh(0.9));
total_time = 5;
dt = 1e-3;
max_it = round(total_time / dt);
% max_it = 10;
GridSize = 128;
init_file = sprintf("%s/initial_phi_%d_r1_%d_r2_%d_space_%d.csv",indir,GridSize, r1,r2,space);
phi0 = readmatrix(init_file);
FileName = sprintf("%s/MG_MATLAB_%d_dt_%.2e_Nx_%d_r1_%d_r2_%d_space_%d",outdir,max_it,dt, GridSize, r1,r2,space);
[final_phi,mass_t,E_t] = CahnHilliard_NMG(phi0,...
                                    t_iter = max_it,...
                                    dt = dt,...
                                    solver_iter = 1e4,...
                                    tol = 1e-6,...
                                    m = 8,...
                                    boundary = "neumann",...
                                    c_relax = 2,...
                                    write_phi = true,...
                                    FileName = FileName);

writematrix(final_phi,sprintf('%s_final_phi.csv', FileName));
writematrix(mass_t,sprintf('%s_mass.csv', FileName));
writematrix(E_t,sprintf('%s_energy.csv', FileName));

fig = figure('visible', 'off');
image(final_phi,'CDataMapping','scaled'); colorbar; axis square;
set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
clim([-1, 1]);
colormap(redblue(100));
saveas(gca,sprintf('%s_final_phi.png', FileName))
% end

function c = redblue(m)
    %   Adam Auton, 9th October 2009
    
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    
    c = [r g b]; 
end    

