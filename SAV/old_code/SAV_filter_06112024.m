% close all;
    a = 0.5;
% Space parameters
    Lx = 1; para.Lx=Lx; Ly=Lx; para.Ly=Ly;
    Nx = 2^7; para.Nx=Nx; Ny=Nx; para.Ny=Nx;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = hx*(0:Nx-1);           
    y  = hy*(0:Ny-1);
    [xx,yy] = meshgrid(x,y); 

    k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);

% Sarah's para
    h=hx;
    h2=h^2;
    M = 8;
    gam = M * h / (2 * sqrt(2) * atanh(0.9)); % gam is our epsilon
    % max_it=15000;
    % para.T=max_it*dt; T=para.T;

% Interface and energy parameters
    epsilon=gam/2^-0; para.epsilon=epsilon;
    C0=0; para.C0=C0;
    M=epsilon^2; para.M=M;
    para.a = a;

% Relaxation and Stabilization
    para.Beta=0; para.gamma0=0;

% Time parameters
    % para.dt   = 1e-1*h2;dt = para.dt;
    % para.dt = 1e-5; dt = para.dt;
    para.dt = 6.25e-6; dt = para.dt;
    para.dtout = 50;  dtout = para.dtout;

    l=1;
    % l=10;
    % l=20;
    % l=40;
    % l=100;
    % l=300;
    % l = 2000;
    % l = 4000;
    % l = 6000;
    para.T=l*dtout*dt; T=para.T;

% Initial data - by function
    T0 = 0; Psi0 = []; % Initialization

    % Tanh function
    R0 = 0.12;
    R = sqrt((xx-0.5).^2 + (yy-0.5).^2);
    eps_c = gam; 
    delta = eps_c*sqrt(2);
    psi0 = 0.5 * (1 - tanh((R - R0)/(2*delta)));
    phi0 = 2*psi0-1;    % psi0=(phi0+1)/2;
    para.phi0 = phi0;

    % phi0 = 0.05*sin(xx).*sin(yy); para.phi0=phi0;
    % phi0 = 0.05*sin(x)'.*sin(y); para.phi0=phi0;
    % r0 = r0_fun(phi0); para.r0=r0;

    % C_cond=10;
    % C_ref =10;
    % psi0=psi0 * C_ref/C_cond;
 
% Initial data - by importing
    filename = "2024-04-23_SAV_output_R0_0.15_eps_0.030019_T_1.8311" + ".mat";
        T0_str = regexp(filename, 'T_([\d.]+(?=\.mat))', 'tokens');
        if ~isempty(T0_str)
            T0 = str2double(T0_str{1}{1});
        else
            error('The value of T0 could not be found in the filename.');
        end
    data = load(filename); 
    Psi0 = data.Psi; Phi0 = 2*Psi0-1;
    phi0 = Phi0(end-Nx+1:end,:);
    para.phi0 = phi0;

    % phi0 = readmatrix('128_100_phi.csv'); % corresponds to C simulation
    % phi0 = load('phi_32_1.0e-5_SD_initial.txt'); para.phi0=phi0;
    % phi0 = load('phi_32_initial.txt'); para.phi0=phi0;
    % phi0 = load('phi_128_initial.txt'); para.phi0=phi0;

% Checking the imaginary part of Fourier coefficients of initial data
    % norm(imag(fft2_filtered(phi0)))
    % norm(imag(fft(phi0(1,:))))
    % norm(imag(fft(phi0(end,:))))
    % norm(imag(fft(phi0(:,1))))
    % norm(imag(fft(phi0(:,end))))

% Filtering out the initial data
    % phi0=ifft2(real(fft2(phi0)));

    % figure;
    % contourf(x,y,phi0); colorbar; axis square

    % figure;
    % err=phi0-ifft2(real(fft2(phi0)));
    % contourf(x,y,err); colorbar; axis square

% Plot initial
    % figure
    % contourf(x,y,psi0); colorbar; axis square
    %     % contourf(x,y,phi0); colorbar; axis square
    %     % contour(x,y,psi0,[0.3,0.3],'k',LineWidth=1.5); axis square
    % set(gca,'FontSize',16);title('t = 0'); xlabel('x'); ylabel('y');

%% CN Solver (Phi collected every dtout results)
    tic;[phi,~,E,~,Phi,O]=CH2d_SAV_CN(para);toc;
    psi=(phi+1)/2; Psi=(Phi+1)/2; 

    % figure; contourf(x,y,psi); colorbar; axis square
    % set(gca,'FontSize',16);title('t = 0'); xlabel('x'); ylabel('y');

%% Post-processing
% plotting for critical radius
    fig_out = 1; % number of figure outputs after initial
    t_Psi   = 0; every = floor(l/fig_out);
    for i = 0:every:l
        T = T0 + i*dtout*dt;
        t_Psi(i/every+1) = T;
        t = Nx*i+1 : Nx*(i+1);
        psi = Psi(t,:);

        figure;
        contourf(x,y,psi); colorbar; axis square;
        set(gca,'FontSize',16);title(['t = ',num2str(T)]); xlabel('x'); ylabel('y');
        filename = sprintf('Contourf_R0_%.2d_N_%d_eps_%.3f_T_%.4f_dt_%.2e.png',R0,Nx,epsilon,T,dt);saveas(gcf, filename);
        close gcf
    end
    if exist('data','var')
        t_Psi = [data.t_Psi t_Psi];
    end


    tt=0;rr=0;everyR=1;
    for i = 0:everyR:l
        T=i*dtout*dt;
        t = Nx*i+1 : Nx*(i+1);
        psi = Psi(t,:);

        % Find and plot the 0.5 contour
        figure('visible', 'off');
        [~, h] = contour(x, y, psi, [0.5 0.5]);
        % Extract the x and y coordinates of the 0.5 contour
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

    if exist('data','var')
        tt = [data.tt tt];
        rr = [data.rr rr];
    end

    % figure(1);
    % % fig_R = figure;
    % 
    % plot(tt,rr,'LineWidth',1.5); hold on
    % xlabel('t');ylabel('R');
    % title('Radius')
    % set(gca,'FontSize',16);
    % filename = sprintf('Radius_evolution_R0_%.2d_N_%d_eps_%.3d_a_%.2d_T_%.2d_dt_%.2e.fig',R0,Nx,epsilon,a,para.T,dt);
    % saveas(gcf, filename);
    % filename = sprintf('Radius_evolution_R0_%.2d_N_%d_eps_%.3d_a_%.2d_T_%.2d_dt_%.2e.png',R0,Nx,epsilon,a,para.T,dt);
    % saveas(gcf, filename);
% 
% %% Plot the modified energy
    dt = para.dt;T = T0 + para.T; tt_E = T0:dt:T;

    if exist('data','var')
        tt_E = [data.tt_E tt_E];
        E    = [data.E E];
    end
% 
%     figure(2);
%     % fig_E = figure;
% 
%     plot(tt_E,epsilon*E,'LineWidth',1.5);set(gca,'FontSize',16); hold on
%     title(['Modified energy, \epsilon=',num2str(epsilon),', dt=',num2str(dt)])
%     xlabel('t');ylabel('E');
%     % hold on;
%     filename = sprintf('Modified_energy_R0_%.2d_N_%d_eps_%.3d_T_%.2d_dt_%.2e.png',R0,Nx,epsilon,para.T,dt);
%     saveas(gcf, filename);
% 
%     % close all;
% 
    current_date = datestr(now, 'yyyy-mm-dd');
    filename = [current_date,'_SAV_output',...
                            '_R0_', num2str(R0),...
                            '_eps_', num2str(epsilon),...
                            '_T_', num2str(T),...
                            '.mat'];
    Psi = [Psi0; Psi];
    save(filename,'Psi','t_Psi','tt','rr','tt_E','E','dt','Nx','Ny','epsilon');


%% Plot the error
    % figure
    %     plot((0:dt:T),D, ...
    %                                 ':', ...
    %         'LineWidth',1.5);set(gca,'FontSize',16);
    %     title('Relative error of r from E1')
    %     xlabel('t');ylabel('(r-sqrt(E1))/sqrt(E1)');
%% Plot the Imaginary magnitude of fft2(phi)
    % figure
    % plot((0:dt:T),O, ...
    %                             ':', ...
    %     'LineWidth',1.5);set(gca,'FontSize',16);
    % title('Imaginary magnitude of fft2(phi)')
    % xlabel('t');ylabel('|imag(fft(phi))|');

%% Video
 % v = VideoWriter('phi0_128_dt_6.25e-6_heatmap.mp4','MPEG-4'); 
 %    % v = VideoWriter('C_Cond_10.mp4','MPEG-4'); % Create a VideoWriter object.
 %    v.FrameRate = 10; % Set to 10 frames per second
 %    open(v); % Open the video file for writing.
 % 
 %    psi0=(phi0+1)/2;
 %    fig = figure;
 % 
 %    heatmap(phi0, 'CellLabelColor','none', 'GridVisible','off');
 %    % contourf(x,y,psi0); axis square; colorbar; % Your plotting code here
 % 
 %    % clim([0 0.35]);
 %    % mesh(xx,yy,phi0,'FaceColor','interp','EdgeColor','interp');colorbar; zlim([-1,1]);% Your plotting code here
 %    set(gca,'FontSize',16);title('t = 0');
 %    frame = getframe(fig); % Capture the current figure window as a frame
 %    writeVideo(v, frame); % Add the frame to the video
 % 
 %    dtout =1;
 %    for i=1:l
 %        T=i*dtout*dt;
 %        t = Nx*i+1 : Nx*(i+1);
 %        psi = Psi(t,:);
 %        fig = figure;
 %        % fig = figure('Visible', 'off');
 % 
 %        % contourf(x,y,psi); colorbar; axis square
 % 
 %        phi=2*psi-1;
 %        heatmap(psi, 'CellLabelColor','none', 'GridVisible','off');
 % 
 %        set(gca,'FontSize',16);title(['t = ',num2str(T)]); xlabel('x'); ylabel('y');
 %        frame = getframe(fig); % Capture the current figure window as a frame
 %        writeVideo(v, frame); % Add the frame to the video
 %        if mod(i,20)==0
 %            close all;
 %        end
 %    end
 %    close(v);

%% Crank Nicolson scheme
function [phi,r,E,D,Phi,O]=CH2d_SAV_CN(para)
% Solve 2D Cahn Hilliard equaiton
% \phi_t = \Delta \mu
% \mu = \Delta\phi + (\phi^3 - \phi)/epsilon^2
global dt epsilon k2 k4 C0 hx hy Lx Ly gamma0 Beta M a

%% Parameters
    T  = para.T;
    dt = para.dt;
    t  = 0;
    Nt = round(T/dt);
    E  = zeros(1,Nt+1);  
    
    Lx = para.Lx;
    Ly = para.Ly;
    Nx= para.Nx;
    Ny= para.Ny;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x = hx*(0:Nx-1);
    y = hy*(0:Ny-1);
    k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/Lx);
    k_y = 1i*[0:Ny/2 -Ny/2+1:-1]*(2*pi/Ly);
    % [kx,ky] = meshgrid(k_x,k_y);
    k_xx = k_x.^2;
    k_yy = k_y.^2;
    [kxx,kyy] = meshgrid(k_xx,k_yy);
    k2 = kxx + kyy;
    % k2 = k_xx+k_yy;
    k4 = k2.^2;
    
    epsilon = para.epsilon;
    C0 = para.C0;
    M = para.M;
    gamma0 = para.gamma0;
    Beta = para.Beta;
    a = para.a;

    phi0 = para.phi0;
    % phi  = phi0;
    % r0 = para.r0;
    r0 = r0_fun(phi0);
    % r  = r0;

    dtout=para.dtout;
    Phi(1:Nx,1:Ny)=phi0;

% Prepare [phi1 r1], the 1st step by FO and CN_pc1st
    % para1 = para; para1.T=dt;
    % [phi1,r1,E,D]  = CH2d_SAV_FO(para1);
    % [phi1,r1,E] = CH1d_SAV_CN_pc1st(para1,phi1);

% Initial modified energy
    E(1)=hx*hy*sum(sum(1/2*phi0.*(-Lap(phi0))))+r0^2;

% Error
    D(1) = (r0-sqrt(hx*hy*sum(sum(f(phi0)))))/sqrt(hx*hy*sum(sum(f(phi0))));
    O(1) = norm(imag(fft2(phi0)));
%% Time iteration
    for nt = 1:Nt
    % for nt = 2:Nt
        t = t+dt;
    
        phi_bar = A_inv_CN(phi0+dt/2*M*Lap(df(phi0)));
        % phi_bar = A_inv_CN(phi1+dt/2*Lap(df(phi1)));
        % phi_bar = 1/2*(3*phi1-phi0);
        
        % Step 1
        b = b_fun(phi_bar);
        
        g = g_fun_CN(phi0,r0,b);
        % g = g_fun_CN(phi1,r1,b);
        
        AiLb = A_inv_CN(M*Lap(b));
        Aig  = A_inv_CN(g);    
        
        gamma = -fft2_filtered(b.*AiLb);
        gamma = gamma(1,1)*hx*hy;
        
        % Step 2      
        bphi = fft2_filtered(b.*Aig);
        bphi = bphi(1,1)*hx*hy/(1+dt/4*gamma);
        
        % Step 3
        phi = dt/4*bphi.*AiLb + Aig;
        r   = r_fun(phi,phi0,r0,b);
        % r   = r_fun(phi,phi1,r1,b);
    
%% update phi0, phi1, r1
        phi0 = phi;
        r0   = r;
        % phi0 = phi1;
        % phi1 = phi;
        % r1  = r;

        if mod(nt,dtout)==0
            ll = nt/dtout;
            tt = Nx*ll+1 : Nx*(ll+1);
            Phi(tt,1:Ny)=phi;
        end

%% Modified energy
        E(1+nt)=hx*hy*sum(sum(1/2*phi.*(-Lap(phi))))+r^2;
% Error
        D(1+nt)=(r-sqrt(hx*hy*sum(sum(f(phi)))))/sqrt(hx*hy*sum(sum(f(phi))));
        O(1+nt) = norm(imag(fft2(phi)));
    end

end

%% General Function Library
function r0 = r0_fun(phi0)
global hx hy C0
    r0 = sqrt(hx*hy*sum(sum(f(phi0))) + C0);
end

function r = r_fun(phi,phi0,r0,b)
global hx hy C0 Beta dt
    bphi0 = fft2_filtered(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);
    bphi  = fft2_filtered(b.*phi);
    bphi  = hx*hy*bphi(1,1);

    E1 = fft2_filtered(f(phi0));
    r = r0 + 1/2*(bphi-bphi0)-Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0)); %r explicit (r multiplied)
end

function b = b_fun(phi)
global C0 hx hy
    E1 = fft2_filtered(f(phi));
    b = df(phi)./sqrt(E1(1,1)*hx*hy+C0);
end

function Lphi=Lap(phi)
global k2
    Lphi=real(ifft2((k2.*fft2_filtered(phi))));
end

% function Dphi = Diff(phi)
% global k
%     Dphi=real(ifft((k.*fft(phi))));
% end

function fphi = f(phi)
global epsilon gamma0 a
    % fphi = (phi.^2-1-gamma0).^2/(4*epsilon^2);
    % a=0.5;
    psi=(phi+1)/2;
    fphi = 4*psi.^2.*(6*a-4*(1+a)*psi+3*psi.^2)/(3*epsilon^2);
    % fphi = (4*psi.^2.*(6*a-4*(1+a)*psi+3*psi.^2)-(4*a-2))/(3*epsilon^2);
end

function dfphi = df(phi)
global epsilon gamma0 a
    % dfphi = (phi.^3-(1+gamma0)*phi)/epsilon^2;

    % a=0.5;
    dpsidphi=1/2;
    psi=(phi+1)/2;
    dfphi = 16*(1-psi).*(a-psi).*psi/epsilon^2 * dpsidphi;
end

%% First Order Function Library
function g = g_fun_FO(phi0,r0,b)
global dt hx hy C0 Beta M
    bphi0 = fft2_filtered(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);

    %r explicit (r multiplied)
    E1 = fft2_filtered(f(phi0));
    g = phi0 + dt.*M*Lap(b)*(r0-1/2*bphi0-Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0))); 
end

function Ai = A_inv_FO(phi)
global dt k2 k4 gamma0 epsilon M
    Ai = real(ifft2(fft2_filtered(phi)./(1+dt*M*k4-dt*gamma0/epsilon^2*M*k2)));
end
%% Crank Nicolson Function Library
function g = g_fun_CN(phi0,r0,b)
global dt hx hy gamma0 epsilon Beta C0 M
    bphi0 = fft2_filtered(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);

    %r explicit (r multiplied)
    E1 = fft2_filtered(f(phi0));
    g = phi0 - dt/2*M*Lap(Lap(phi0)) + dt/2*gamma0/epsilon^2*M*Lap(phi0) ...
    + dt*M.*Lap(b)*(r0-1/4*bphi0-1/2*Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0)));
end

function Ai = A_inv_CN(phi)
global dt k2 k4 gamma0 epsilon M
    Ai = real(ifft2(fft2_filtered(phi)./(1+dt/2*M*k4-dt/2*gamma0/epsilon^2*M*k2)));
end

function y=fft2_filtered(x)
    y=real(fft2(x));
end
