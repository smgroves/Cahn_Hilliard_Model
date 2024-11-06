
a = 0.5;
nt0 = 0; para.nt0 = nt0;

% Arrays to store Nx and elapsed times
Ns = 2.^(10);
dt_factors = [0.25 0.5 1 2 4]; 
% dt_factors = [5];
times = zeros(length(Ns), length(dt_factors));

% For Error Curves Plotting
% Ns = 2.^(9);
% dt_factors = [0.25 0.5 1 2 4]; 
% Ds = cell(length(Ns), length(dt_factors)); % To store error data
% Ts = cell(length(Ns), length(dt_factors)); % To store time vectors for error

% Original dt and T_total
original_dt = 2.5e-5;
dtout = 1e1; 
l_original = 120; 
T_total = l_original * dtout * original_dt; % Total simulation time remains constant 0.03

% Loop over Nx from 2^6 to 2^11
for idx = 1:length(Ns)   
 % Space parameters
    Lx = 2; para.Lx=Lx; Ly=Lx; para.Ly=Ly;
    Nx = Ns(idx); para.Nx=Nx; Ny=Nx; para.Ny=Nx;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = hx*(0:Nx/2-1);           
    y  = hy*(0:Ny/2-1);
    [xx,yy] = meshgrid(x,y); 

% Sarah's para
    h=hx;
    h2=h^2;
    M = 8;
    gam = M * h / (2 * sqrt(2) * atanh(0.9)); % gam is our epsilon

% Interface and energy parameters
    epsilon=1*gam; para.epsilon=epsilon;
    C0=0; para.C0=C0;
    M=epsilon^2; para.M=M;
    para.a = a;

% Relaxation and Stabilization
    para.Beta=0; Beta=para.Beta;
    para.gamma0=0;

% Time parameters
    para.dtout = dtout;
 
% Initial data - by importing
   Nphi = Nx / 2;
   filename = sprintf('initial_phi_%d_mean_0_sd_0.2.csv', Nphi);
   phi0 = readmatrix(filename);

% For each dt scaling factor
    for dt_idx = 1:length(dt_factors)
        dt_factor = dt_factors(dt_idx);
        dt = dt_factor * original_dt; para.dt = dt;

        % Adjust l to keep T constant
        % l_adjusted = round(T_total / (dtout * dt));
        l_adjusted = 480;
        para.T = l_adjusted * dtout * dt; T = para.T;

        para.phi0 = phi0;
        phi0 = ext(phi0);   
        para.phi0 = phi0;
    
    %% CN Solver (Phi collected every dtout results)
        tic;
        [phi,~,E,D,Phi,O]=CH2d_SAV_CN(para);
        time_elapsed = toc;
    
        % Store the elapsed time
        times(idx, dt_idx) = time_elapsed;

        % % Store the error data and time vector
        % t_vec = 0:dt:para.T;
        % Ds{idx, dt_idx} = D;
        % Ts{idx, dt_idx} = t_vec;
    
        phi0 = extback(phi0); 
        phi = extback(phi);
    
        psi=(phi+1)/2; Psi=(Phi+1)/2; 
        
        Lx = Lx/2;Ly=Lx;
        Nx = Nphi;Ny=Nx;
        
        hx = Lx/Nx;
        hy = Ly/Ny;
        x  = hx*(0:Nx-1);           
        y  = hy*(0:Ny-1);
    
        current_date = datestr(now, 'yyyy-mm-dd');
        filename = [current_date,'_SAV_',...
                                    'SD',...
                                    '_a_', num2str(a),...
                                    '_eps_', num2str(epsilon),...
                                    '_N_', num2str(Nx),...
                                    '_dt_', num2str(dt),...
                                    '_T_', num2str(T),...
                                    '_Beta_', num2str(Beta),...
                                    '_C0_', num2str(C0),...
                                    '_nt0_', num2str(nt0)];
        filename_mat = [filename, '.mat'];
        save(filename_mat,'Psi','E','dt','Nx','Ny','epsilon');
        
        % Display progress
        fprintf('Nx = %d , dt = %d completed in %f seconds.\n', Nphi, dt, time_elapsed);
    end
end

% % Display the times matrix
% disp('Elapsed times for each Nx and dt scaling factor:');
% disp('Rows correspond to Nx values, columns to dt factors');
% disp('Nx values:');
% disp(Ns/2');
% disp('dt scaling factors:');
% disp(dt_factors);
% disp('Elapsed times (in seconds):');
% disp(times);

% figure; % Plot CPU time vs. Nx for each dt scaling factor
% for dt_idx = 1:length(dt_factors)
%     dt_factor = dt_factors(dt_idx);
%     loglog(Ns/2, times(:, dt_idx), '-o', 'LineWidth', 2, 'DisplayName', sprintf('dt factor = %.2f', dt_factor));
%     hold on;
% end

% xlabel('Nx (Spatial Resolution)');
% ylabel('CPU Time (seconds)');
% title('CPU Time vs. Nx for Different dt Scaling Factors');
% legend('Location', 'NorthWest');
% grid on;

% % Adjust the x-ticks to show powers of 2
% set(gca, 'XTick', Ns);
% set(gca, 'XTickLabel', arrayfun(@num2str, Ns, 'UniformOutput', false));

% filename = sprintf('cpu_time_vs_nx_varying_dt.png');
% saveas(gcf, filename);


%% Plot Error Curves
% figure;hold on;
% % Define colors for each curve
% num_curves = length(Ns) * length(dt_factors);
% colors = lines(num_curves); % Generate distinguishable colors
% legend_entries = cell(num_curves, 1);
% legend_idx = 1;
% 
% for idx = 1:length(Ns)
%     Nx = Ns(idx);
%     for dt_idx = 1:length(dt_factors)
%         dt_factor = dt_factors(dt_idx);
%         dt = dt_factor * original_dt;
% 
%         % Retrieve D and t_vec
%         D = Ds{idx, dt_idx};
%         t_vec = Ts{idx, dt_idx};
% 
%         % Plot the error curve
%         plot(t_vec, D, 'LineWidth', 1.5, 'Color', colors(legend_idx, :), ...
%             'DisplayName', sprintf('Nx=%d, dt=%.4g', Nx/2, dt));
% 
%         legend_idx = legend_idx + 1;
%     end
% end
% 
% set(gca, 'FontSize', 16);
% title('Relative Error of r from E1');
% xlabel('t');
% ylabel('(r - sqrt(E1)) / sqrt(E1)');
% legend('show', 'Location', 'BestOutside');
% grid on;
% hold off;
% set(gcf, 'Position', [100, 100, 800, 600]); 
% 
% filename = sprintf('error_curves_Nx_%d.png', Nx/2);
% saveas(gcf, filename);


%% Crank Nicolson scheme
function [phi,r,E,D,Phi,O]=CH2d_SAV_CN(para)
% Solve 2D Cahn Hilliard equaiton
% \phi_t = \Delta \mu
% \mu = \Delta\phi + (\phi^3 - \phi)/epsilon^2
global dt epsilon k2 k4 C0 hx hy Lx Ly gamma0 Beta M a nt nt0

%% Parameters
    nt0 = para.nt0;

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
        
        gamma = -fft2(b.*AiLb);
        gamma = gamma(1,1)*hx*hy;
        
        % Step 2      
        bphi = fft2(b.*Aig);
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
    bphi0 = fft2(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);
    bphi  = fft2(b.*phi);
    bphi  = hx*hy*bphi(1,1);

    E1 = fft2(f(phi0));
    r = r0 + 1/2*(bphi-bphi0)-Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0)); %r explicit (r multiplied)
end

function b = b_fun(phi)
global C0 hx hy
    E1 = fft2(f(phi));
    b = df(phi)./sqrt(E1(1,1)*hx*hy+C0);
end

function Lphi=Lap(phi)
global k2
    Lphi = real(ifft2((k2.*fft2_filtered(phi))));
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
    bphi0 = fft2(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);

    %r explicit (r multiplied)
    E1 = fft2(f(phi0));
    g = phi0 + dt.*M*Lap(b)*(r0-1/2*bphi0-Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0))); 
end

function Ai = A_inv_FO(phi)
global dt k2 k4 gamma0 epsilon M
    Ai = real(ifft2(fft2_filtered(phi)./(1+dt*M*k4-dt*gamma0/epsilon^2*M*k2)));
end
%% Crank Nicolson Function Library
function g = g_fun_CN(phi0,r0,b)
global dt hx hy gamma0 epsilon Beta C0 M
    bphi0 = fft2(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);

    %r explicit (r multiplied)
    E1 = fft2(f(phi0));
    g = phi0 - dt/2*M*Lap(Lap(phi0)) + dt/2*gamma0/epsilon^2*M*Lap(phi0) ...
    + dt*M.*Lap(b)*(r0-1/4*bphi0-1/2*Beta*dt*r0*(r0-sqrt(E1(1,1)*hx*hy+C0)));
end

function Ai = A_inv_CN(phi)
global dt k2 k4 gamma0 epsilon M
    Ai = real(ifft2(fft2_filtered(phi)./(1+dt/2*M*k4-dt/2*gamma0/epsilon^2*M*k2)));
end

function y = fft2_filtered(x)
global nt nt0
    if nt>=nt0
        % y=real(fft2(x));
        y=fft2(x);
    else
        y=fft2(x);
    end
end

function x_ext = ext(x)
    N = length(x);
    x_ext = zeros(2*N);
    x_ext(1:N, 1:N)         = x;
    x_ext(1:N, N+1:2*N)     = x(:       , end:-1:1);
    x_ext(N+1:2*N, 1:N)     = x(end:-1:1,    :    );
    x_ext(N+1:2*N, N+1:2*N) = x(end:-1:1, end:-1:1);
end

function x_extback = extback(x)
    N = length(x);
    x_extback=x(1:N/2,1:N/2);
end