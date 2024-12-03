indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/";
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function";
m = 8;
GridSize = 128;
% dt = .9e-2/(GridSize^2);
dt = 5.4e-7;
n_relax = 4;
h = 1/GridSize;
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
% epsilon = m/(2 * sqrt(2) * atanh(0.9));
gamma = epsilon^2/h^2;
D = GridSize^2;
% total_time = (.1/64^2)*500;
total_time = dt*2000;
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
image(final_phi,'CDataMapping','scaled'); colorbar; axis square;
set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
clim([-1, 1]);
colormap(redblue(100));
saveas(gca,sprintf('%s_final_phi.png', FileName))

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


% /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "run_fd_droplets();quit;"
