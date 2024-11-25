indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/";
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output";

for GridSize = [128]
    for dt = [ 1e-2/(GridSize^2), 1e-1/(GridSize^2)]
        for n_relax = [2, 4, 8]
            m = 8;
            h = 1/GridSize;
            epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
            % epsilon = m/(2 * sqrt(2) * atanh(0.9));
            total_time = (.1/64^2)*500;
            max_it = round(total_time / dt);
            % max_it = 4000;
            % GridSize = 128;
            FileName = sprintf("%s/MG_%d_dt_%.9e_Nx_%d_n_relax_%d_eps_0.015009369912862116",outdir,max_it,dt, GridSize,n_relax);
            FileName = regexprep(FileName, 'e(-?)0+', 'e$1')
            phi = readmatrix(sprintf("%s_final_phi.csv",FileName));


            fig = figure('visible', 'off');
            image(phi,'CDataMapping','scaled'); colorbar; axis square;
            set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
            clim([-1, 1]);
            colormap(redblue(100));
            saveas(gca,sprintf('%s_final_phi.png', FileName))
        end
    end
end

% outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/checkerboard_IC/output";
% GridSize = 128;
% dt =  1/(GridSize^2);
% m = 8;
% h = 1/GridSize;
% epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9));
% total_time = (1/64^2)*5000;
% max_it = round(total_time / dt);
% FileName = sprintf("%s/MG_%d_dt_%.9e_Nx_%d_eps_0.015009369912862116_grid_8",outdir,max_it,dt, GridSize);
% FileName = regexprep(FileName, 'e(-?)0+', 'e$1')
% phi = readmatrix(sprintf("%s_final_phi.csv",FileName));
% fig = figure('visible', 'off');
% image(phi,'CDataMapping','scaled'); colorbar; axis square;
% set(gca,'FontSize',16);title(['t = ',num2str(total_time)]); xlabel('x'); ylabel('y');
% clim([-1, 1]);
% colormap(redblue(100));
% saveas(gca,sprintf('%s_final_phi.png', FileName))


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


