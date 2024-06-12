dt_name = "6.25e-6";
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet";
total_time=0.0625;
dt = str2double(dt_name);
timesteps=total_time/dt;
name=sprintf('phi_128_%s_1.0e-6__dt_%s',string(timesteps), dt_name);
dtout=10;
Nx = 2^7;Ny=2^7;
hx = 1/Nx;
hy = 1/Ny;
x  = hx*(0:Nx-1);           
y  = hy*(0:Ny-1);
[xx,yy] = meshgrid(x,y); 
T0=0;
tt=0;rr=0;everyR=10;
phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
% phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
phidims = size(phi);
phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
phidims(1) = phidims(2); %Determine size of square grid
phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension

for i = 0:everyR:timesteps/dtout
    T=i*dtout*dt;
    % t = Nx*i+1 : Nx*(i+1);
    % psi = Psi(t,:);
    phi_tmp = phi(:,:,i+1);
    % Find and plot the 0.5 contour
    figure('visible', 'off');
    [~, h] = contour(x, y, phi_tmp, [0 0]);
    % Extract the x and y coordinates of the 0 contour
    contour_data = h.ContourMatrix;
    x_contour = contour_data(1, 2:end); 
    y_contour = contour_data(2, 2:end);
    distances = sqrt((x_contour-0.5).^2 + (y_contour-0.5).^2);
    radius = mean(distances);
    close
    tt(i/everyR+1) = T0 + T;
    rr(i/everyR+1) = radius;
end
% clear h; 

%% Plot the figure
    fig = figure;
        plot(tt,rr, ...
                                    ':', ...
            'LineWidth',1.5);set(gca,'FontSize',16);
        title(sprintf('Radius at 0.5 level set %s',name))
        xlabel('t');ylabel('R');
        filename = sprintf('%s/radius_0.5_nx_128_nt_%s_1.0e-6__dt_%s.png',indir,string(timesteps),dt_name);
        saveas(fig, filename);