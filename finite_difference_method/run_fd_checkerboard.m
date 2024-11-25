
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/checkerboard_IC/IC/";
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/checkerboard_IC/";
grid = 8
m = 8;
GridSize = 128
h = 1/GridSize;
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
% epsilon = m/(2 * sqrt(2) * atanh(0.9));
gamma = epsilon^2/h^2;
D = GridSize^2;
total_time = (1/64^2)*500;
dt_orig = .1 / (GridSize^2);
for factor = [.1]
% total_time = 0.01;
    dt = dt_orig*factor
    max_it = round(total_time / dt);
    FileName = sprintf("%s/FD_%d_dt_%.2e_Nx_%d_gam_%.2e_D_%g_grid_%d",outdir,max_it,dt, GridSize,gamma,D, grid)

    [final_phi,mass_t,E_t] = spinodal_decomp(D,gamma, ...
                    "dt",dt,...
                    "GridSize",GridSize,...
                    "CaptureMode","standard",...
                    "NumIterations",max_it,...
                    "FrameSpacing",100/factor,...
                    'ImgStyle','true',...
                    "InputMatrix",sprintf("%s/initial_phi_%d_grid_%d.csv",indir,GridSize, grid),...
                    "FileName",FileName,...
                    "InputType",'phi',...
                    "ConstantColorbar", true,...
                    "write_residual",false,...
                    "write_phi",false);
            % save([FileName, '.mat'],'phi_t');
    writematrix(final_phi,sprintf('%s_final_phi.csv', FileName));
    writematrix(mass_t,sprintf('%s_mass.csv', FileName));
    writematrix(E_t,sprintf('%s_energy.csv', FileName));

    fig = figure('visible', 'off');
    image(final_phi,'CDataMapping','scaled'); colorbar; axis square;
    set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
    clim([-1, 1]);
    colormap(redblue(100));
    saveas(gca,sprintf('%s_final_phi.png', FileName))
end

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
