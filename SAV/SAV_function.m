% This code has been adapted by SM Groves on November 25, 2024, to be called in external scripts. 
function [final_psi, E] = run_CH_solver(Psi0, m, t_iter, dt, options)
    arguments
    Psi0 = "";
    m = 8;
    t_iter = 100;
    dt = 1e-5;
    options.Lx double = 1.0;
    options.Ly double = 1.0;
    options.dt_out int = 10;
    end
% close all;
    a = 0.5; %Needs to be removed because it is always 0.5 now
% Space parameters
    para.Lx=options.Lx; para.Ly=options.Ly;
    [Nx,Ny] = size(phi0); %Define number of grid points in x and y
    para.Nx=Nx; para.Ny=Nx;
    hx = options.Lx/Nx;
    hy = options.Ly/Ny;
    x  = hx*(0:Nx-1);           
    y  = hy*(0:Ny-1);
    [xx,yy] = meshgrid(x,y); 

    k_x = 1i*[0:Nx/2 -Nx/2+1:-1]*(2*pi/options.Lx);

% Sarah's para
    h2=hx*hy;
    epsilon = m * h / (2 * sqrt(2) * atanh(0.9)); % gam is our epsilon

% Interface and energy parameters
    para.epsilon=epsilon;
    C0=0; para.C0=C0;
    M=epsilon^2; para.M=M;
    para.a = a;

% Relaxation and Stabilization
    para.Beta=0; para.gamma0=0;

% Time parameters
    % para.dt   = 1e-1*h2;dt = para.dt;
    % para.dt = 1e-5; dt = para.dt;
    para.dt = dt;
    para.dtout = options.dt_out;  

    para.T=t_iter*dt; T=para.T;

    Phi0 = 2*Psi0-1;
    phi0 = Phi0(end-Nx+1:end,:);
    para.phi0 = phi0;


%% CN Solver (Phi collected every dtout results)
    tic;[phi,~,E,~,Phi,O]=CH2d_SAV_CN(para);toc;
    psi=(phi+1)/2; Psi=(Phi+1)/2; 

%% Post-processing

% %% Plot the modified energy
    % dt = para.dt;T = T0 + para.T; tt_E = T0:dt:T;

    % if exist('data','var')
    %     tt_E = [data.tt_E tt_E];
    %     E    = [data.E E];
    % end

    % figure(2);
    % % fig_E = figure;

    % plot(tt_E,epsilon*E,'LineWidth',1.5);set(gca,'FontSize',16); hold on
    % title(['Modified energy, \epsilon=',num2str(epsilon),', dt=',num2str(dt)])
    % xlabel('t');ylabel('E');
    % % hold on;
    % filename = sprintf('Modified_energy_R0_%.2d_N_%d_eps_%.3d_T_%.2d_dt_%.2e.png',R0,Nx,epsilon,para.T,dt);
    % saveas(gcf, filename);

    % % close all;

    current_date = datestr(now, 'yyyy-mm-dd');
    filename = [current_date,'_SAV_output',...
                            '_R0_', num2str(R0),...
                            '_eps_', num2str(epsilon),...
                            '_T_', num2str(T),...
                            '.mat'];
    Psi = [Psi0; Psi];
    save(filename,'Psi','t_Psi','tt','rr','tt_E','E','dt','Nx','Ny','epsilon');

end


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
