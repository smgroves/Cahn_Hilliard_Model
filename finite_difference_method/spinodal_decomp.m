%source: https://github.com/nsbalbi/Spinodal-Decomposition
% adapted to use same initial conditions as our method 
% and to output energy, mass, and residual

function spinodal_decomp(D,gamma,options)
% SPINODAL_DECOMP Generates and records Cahn-Hilliard spinodal
% decomposition models using Euler's method.
% 
%    SPINODAL_DECOMP() generates a Cahn-Hilliad spindoal decomposition
%    model with sample coefficients. The result is saved as an AVI named
%    'spindoal_decomposition'.
%
%    SPINODAL_DECOMP(D,gamma) generates a Cahn-Hilliad spindoal 
%    decomposition model with diffusion coefficient D and square length of
%    transitional regions gamma. The result is saved as an AVI named 
%    'spindoal_decomposition'. See "Inputs" for details.
%    
%    SPINODAL_DECOMP(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
%    name-value pairs that can be used to change default function
%    parameters. See "Parameters" for details.
%
% Inputs
%    D: double. Diffusion coefficient. Default is 10
%    gamma: double. Sqaure length of transitional regions between domains.
%       Default is 5
% 
% Parameters
%    'dt': double. Time increment between iterations of Euler's method. 
%       (Default = 0.005)
%    'GridSize': int. Edge length for the square grid used for the
%       model. Note generation time increases exponentially with grid size.
%       (Default = 200)
%    'NumIterations': int. Total number of model iterations. Each iteration
%       represents a step of Euler's method. (Default = 10000)
%    'FrameSpacing': int. Number of iterations between captured frames, 
%       only relevant if capture mode is standard. (Default = 10)
%    'CaptureMode': char. Method of video cature. Possible inputs below.
%         'standard' - Constant num of iterations between frames. (Default)
%      'incremental' - Num iterations between frames increases over time
%    'ImgStyle': char. Method of frame generation. Possible inputs below.
%         'binary' - Concentrations are binarized to major species. (Default)
%           'true' - Concentrations are mapped to the colormap.
%    'Colormap': char. Colormap used for visualization. Supports all 
%       default MATLAB colormaps (Default = 'pink')
%    'FileName': char. Name of video (Default = 'spinodal_decomposition')
%       
% Examples
%    Example 1: Generate and record default model.
%       spinodal_decomp();
%    Example 2: Generate and record model with custom constants.
%       spinodal_decomp(20,3);
%    Example 3: Generate and record model with custom constants and capture
%    mode.
%       spinodal_decomp(10,10,'CaptureMode','incremental');
%    Example 4: Generate and record model with custom constants and
%    multiple custom parameters.
%       spinodal_decomp(10,10,...
%                      'CaptureMode','incremental',...
%                      'Colormap','jet',...
%                      'ImgStyle','true',...
%                      'NumIterations',25000);   

arguments
   D double = 10
   gamma double = 5
   options.dt double = 0.005
   options.GridSize double = 200
   options.NumIterations double = 10000
   options.FrameSpacing double = 10
   options.CaptureMode char = 'standard'
   options.ImgStyle char = 'binary'
%    options.Colormap char = 'pink'
   options.FileName char = 'spinodal_decomposition'
   options.InputMatrix char = 'input.csv'
   options.InputType char = 'phi'
   options.ConstantColorbar = true
end

% Random Starting Concentrations (-1 and 1 represent different species)
% u = 2*randi(2,options.GridSize)-3;
u = readmatrix(options.InputMatrix);
% if options.InputType == 'phi'
%     u = (u+1)./2;
% end


writer = VideoWriter(options.FileName,'MPEG-4');
open(writer);
figure(1)
% colormap(options.Colormap)
colormap(redblue(100))
image(u,'CDataMapping','scaled'); colorbar; axis square;
set(gca,'FontSize',16);title('t = 0'); xlabel('x'); ylabel('y');

frame = getframe(1);

writeVideo(writer,frame);

% Variables for incremental capture
count = 1;
frameStep = 1;

for i = 1:options.NumIterations
   curr_t=(i)*options.dt;

   u = iterate(u,D,gamma,options.dt);
   % Incremental video mode
   if (strcmp(options.CaptureMode,'incremental'))
       if (count == frameStep)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');colorbar; axis square;
           else
               uMod = round((u+1)/2); % Binarizes image
               image(uMod,'CDataMapping','scaled');colorbar; axis square;
           end
           set(gca,'FontSize',16);title(['t = ',num2str(curr_t)]); xlabel('x'); ylabel('y');
           frame = getframe(1);
           writeVideo(writer,frame);
           count = 0;
           frameStep = frameStep + 1;
       end
       count = count+1;
   % Standard video mode
   else
       if (mod(i,options.FrameSpacing) == 0)
           if (strcmp(options.ImgStyle,'true'))
               image(u,'CDataMapping','scaled');colorbar; axis square;
           else
               uMod = round((u+1)/2);
               image(uMod,'CDataMapping','scaled');colorbar; axis square;
           end
           if options.ConstantColorbar
               clim([-1, 1]);
           end
           set(gca,'FontSize',16);title(['t = ',num2str(curr_t)]); xlabel('x'); ylabel('y');
           frame = getframe(1);
           writeVideo(writer,frame);
       end
   end
end

close(writer);

fprintf('Done!\n');

% Forward Euler Method iteration of model
function uOut = iterate(u,D,gamma,dt) 
    % Calculates laplacian of concentration field
    uLaplace = laplacian(u);
    % Calculates chemical potentials
    uMu = u.^3 - u - gamma*uLaplace;
    % Laplacian of chemical potentials
    muLaplace = laplacian(uMu);
    % Cahn-Hilliard Equation
    duT = D*muLaplace;
    % Foreward Euler Method
    uOut = u + dt*duT;
    res2 = calc_residual(uOut,u,uMu,dt,D);
    fid = fopen(sprintf('%s_residual.txt',options.FileName), 'a+');
    fprintf(fid, '%f \n', res2);
    fclose(fid);
    E = discrete_energy(uOut,u,uMu,dt,D);
    fid = fopen(sprintf('%s_energy.txt',options.FileName), 'a+');
    fprintf(fid, '%f \n', E);
    fclose(fid);
    % writematrix(uOut,sprintf('%s_phi.csv', options.FileName),'WriteMode','append');
end 

end

function res2 = calc_residual(uNew, uOld, uMu, dt, D)
    s = size(uNew);
    Nx = s(1);
    rr = uMu;
    sor = D*laplacian(rr);
    rr = sor - (uNew - uOld)/dt;
    x = sum(sum(rr.*rr));
    res2 = sqrt(x / (Nx * Nx));
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

function E = discrete_energy(phi,hxy,gridx,gridy,eps2)
    %Local function for calculating the discrete energy across the domain
    f = @(x) 0.25*(x.^2-1).^2; %Define double-well free energy
    a = hxy*sum(sum(f(phi))); %Calculate chemical free energy
    sum_i = 0; %Initialize interfacial free energy in x
    for i = 1:gridx-1
        for j = 1:gridy
            sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
        end
    end
    sum_j = 0; %Initialize interfacial free energy in y
    for i = 1:gridx
        for j = 1:gridy-1
            sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
        end
    end
    E = a + 0.5*eps2*(sum_i+sum_j);
    end