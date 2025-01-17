function [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV(phi0, varargin)
%This function uses the scalar auxiliary variable method to solve the 
%Cahn-Hilliard equation for a specified number of time steps of size dt.
% 
%INPUTS
    %phi0 = Initial field of chemical states in the domain, created by ch_initialization.
%
%NAME-VALUE PAIRS
    %t_iter = Number of time steps simulated. Default = 1e3.
    %dt = Time step. Default = 2.5e-5 characteristic times.
    %solver_iter = Number of solver iterations per time step. Default = 1e4.
    %dt_output = Number of time steps output to phi_t as a multidimensional
    %   array (if less than 1e9 elements) or printed file (if greater than
    %   1e9 elements). Default = nan (output as many timesteps that enable
    %   saving a multidimensional array of less than 1e9 elements).
    %m = Number of mesh points over which the interface exists. Default = 4. **OR SHOULD THIS BE 8**
    %epsilon2 = Squared interfacial transition distance; if specified,
    %   m will be overwritten. Default = nan (do not overwrite m).
    %boundary = Boundary conditions for the simulation:
    %   'periodic' (default) - flux on one domain border equals negative flux on the opposite border.
    %   'neumann' - zero flux on the domain borders.
    %domain = Vector of rightmost and leftmost grid points in x and y.
    %   Format: [xright xleft yright yleft]. Default = [1 0 1 0].
    %printphi = Logical to print phi to a file regardless of whether or 
    %   not it can be saved as a multidimensional array. Default = false.
%
%OUTPUT
    %t_out = Time corresponding to the dt time step outputs.
    %phi_t = Multidimensional array of phi over t_out.
    %delta_mass_t = Vector of mass change over t_out.
    %E_t = Vector of relative energy over t_out.
    
    % Input Parsing and Initialization
    % --------------------------------
    default_t_iter = 1e3;
    default_dt = 2.5e-5;
    default_dt_output = nan;
    default_m = 4;
    default_epsilon2 = nan;
    default_boundary = 'periodic';
    default_domain = [1 0 1 0];
    default_printphi = false;
    default_C0 = 0;
    default_Beta = 0;
    default_gamma0 = 0;
    
    CahnHilliard_SAV_parser = inputParser;

    %Set general criteria for inputs and name-value pairs
    valid_matrix = @(x) ismatrix(x);
    valid_integer = @(x) x-floor(x) == 0;
    valid_pos_num = @(x) isnumeric(x) && (x > 0);
    valid_boundary_type = @(x) strcmpi(x,'periodic') || strcmpi(x,'neumann');
    valid_domain_vector = @(x) length(x) == 4;
    valid_logical = @(x) islogical(x) || x == 1 || x == 0;

    %Set parser options and valid input criteria
    addRequired(CahnHilliard_SAV_parser,'phi0',valid_matrix);
    
    addParameter(CahnHilliard_SAV_parser,'t_iter',default_t_iter,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'dt',default_dt,valid_pos_num);
    addParameter(CahnHilliard_SAV_parser,'dt_output',default_dt_output,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'m',default_m,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'epsilon2',default_epsilon2,valid_pos_num);
    addParameter(CahnHilliard_SAV_parser,'domain',default_domain,valid_domain_vector);
    addParameter(CahnHilliard_SAV_parser,'boundary',default_boundary,valid_boundary_type);
    addParameter(CahnHilliard_SAV_parser,'printphi',default_printphi,valid_logical);
    addParameter(CahnHilliard_SAV_parser,'C0',default_C0,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'Beta',default_Beta,valid_integer);
    addParameter(CahnHilliard_SAV_parser,'gamma0',default_gamma0,valid_integer);

    parse(CahnHilliard_SAV_parser, phi0, varargin{:});
    
    %Extract parsed inputs
    phi0 = CahnHilliard_SAV_parser.Results.phi0;
    t_iter = CahnHilliard_SAV_parser.Results.t_iter;
    dt = CahnHilliard_SAV_parser.Results.dt;
    dt_output = CahnHilliard_SAV_parser.Results.dt_output;
    m = CahnHilliard_SAV_parser.Results.m;
    epsilon2 = CahnHilliard_SAV_parser.Results.epsilon2;
    boundary = CahnHilliard_SAV_parser.Results.boundary;
    xright = CahnHilliard_SAV_parser.Results.domain(1);
    xleft = CahnHilliard_SAV_parser.Results.domain(2);
    yright = CahnHilliard_SAV_parser.Results.domain(3);
    yleft = CahnHilliard_SAV_parser.Results.domain(4);
    printphi = CahnHilliard_SAV_parser.Results.printphi;
    C0 = CahnHilliard_SAV_parser.Results.C0;
    Beta = CahnHilliard_SAV_parser.Results.Beta;
    gamma0 = CahnHilliard_SAV_parser.Results.gamma0;

%% Define and initialize key simulation parameters

[nx,ny] = size(phi0); %Define number of grid points in x and y
%Define maximum possible V-cycle levels as the largest power of two
%shared by both x and y dimensions
i = 1;
while i < min([floor(log2(nx)) floor(log2(ny))]) && mod(nx,2^i) == 0 ...
    && mod(ny,2^i) == 0
    i = i+1;
end
% n_level = i; clear i;
Lx = xright-xleft; Ly = yright-yleft;
hx = Lx/nx; hy = Ly/ny;
h2 = hx*hy; %Define mesh size
x = hx*(0:nx-1); y = hy*(0:ny-1);
k_x = 1i*[0:nx/2 -nx/2+1:-1]*(2*pi/Lx); k_y = 1i*[0:ny/2 -ny/2+1:-1]*(2*pi/Ly);
k_xx = k_x.^2; k_yy = k_y.^2;
[kxx,kyy] = meshgrid(k_xx,k_yy);
k2 = kxx + kyy;
k4 = k2.^2;
if isnan(epsilon2)
    epsilon2 = h2*m^2/(2*sqrt(2)*atanh(0.9))^2; %Define ϵ^2 if not prespecified
else
    m = sqrt((epsilon2*(2*sqrt(2)*atanh(0.9))^2)/h2); %Else overwrite m
end

phi_old = phi0; %Initialize prior chemical state
phi_new = phi0; %Initialize next chemical state
r_old = r0_fun(phi0,hx,hy,C0); %Initialize prior sav state
r_new = r_old; %Initialize next sav state
downsampled = nx*ny*t_iter > 1e9; %Logical index for the need to downsample
optdownsampled = dt_output < t_iter; %Logical index if the user opted to downsample
if ~downsampled && ~optdownsampled %Initialize outputs
    phi_t = zeros(nx,ny,t_iter);
    mass_t = zeros(t_iter,1);
    E_t = zeros(t_iter,1);
    t_out = 0:dt:((t_iter-1)*dt);
    t_spacing = 1;
else %Downsample outputs
    if ~optdownsampled
        t_iter_ds = floor(1e9/nx/ny);
        t_spacing = floor(t_iter/t_iter_ds);
    else
        t_iter_ds = dt_output;
        t_spacing = floor(t_iter/dt_output); %Space out according to dt_output
    end
    if isnan(dt_output) || optdownsampled %If downsampled to be saved
        phi_t = zeros(nx,ny,t_iter_ds);
    end
    mass_t = zeros(t_iter_ds,1);
    E_t = zeros(t_iter_ds,1);
    t_out = 0:t_spacing:(t_iter-t_spacing);
end
    
%% Run SAV solver

if ~downsampled && (isnan(dt_output) || ~optdownsampled || ~printphi) %If output is not specified or does not need to be downsampled
    for i = 1:t_iter
        phi_t(:,:,i) = phi_old;
        mass_t(i) = sum(sum(phi_old))/(h2*nx*ny);
        E_t(i) = ch_discrete_energy(phi_old,h2,nx,ny,epsilon2);
        phi_new = sav_solver(phi_old,phi_new,r_old,r_new,nx,ny,hx,hy,k2,k4, ...
            xright,xleft,yright,yleft,dt,epsilon2, ...
            boundary,C0,Beta,gamma0);
        phi_old = phi_new;
        r_old = r_new;
        if mod(i/t_iter*100,5) == 0
            fprintf('%3.0f percent complete\n',i/t_iter*100)
        end
    end
else %Downsample output or specify output
    if isnan(dt_output)
        fprintf('Downsaving every %4.0f time steps\n',t_spacing)
    else
        fprintf('Saving phi_t to file in the working directory\n')
    end       
    k = 1; %Initialize counter and outputs
    phi_t(:,:,k) = phi_old;
    mass_t(k) = sum(sum(phi_old))/(h2*nx*ny);
    E_t(k) = ch_discrete_energy(phi_old,h2,nx,ny,epsilon2);
    for i = 0:t_spacing:(t_iter-t_spacing)
        for j = 1:t_spacing %Iterate through t_spacing steps
            phi_temp = sav_solver(phi_old,phi_new,r_old,r_new,nx,ny,hx,hy,k2,k4, ...
                xright,xleft,yright,yleft,dt,epsilon2, ...
                boundary,C0,Beta,gamma0);
            phi_old = phi_temp;
            r_old = r_new;
        end
        if isnan(dt_output) && ~printphi %Save to variable
            phi_t(:,:,k) = phi_temp;
        else %Write to file
            writematrix(phi_temp,strcat(pwd,'/phi_t.csv'),'WriteMode','append');
        end
        mass_t(k) = sum(sum(phi_temp))/(h2*nx*ny);
        E_t(k) = ch_discrete_energy(phi_temp,h2,nx,ny,epsilon2);
        k = k+1;
        if mod(i/t_iter*100,5) == 0
            fprintf('%3.0f percent complete\n',i/t_iter*100)
        end
    end
end

%Center mass and normalize energy to t == 0
delta_mass_t = mass_t - mass_t(1);
E_t = E_t/E_t(1);

end
