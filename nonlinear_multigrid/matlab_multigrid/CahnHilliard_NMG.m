function [final_phi, mass_t,E_t] = CahnHilliard_NMG(phi0,varargin)
%This function uses the nonlinear multigrid method to solve the 
%Cahn-Hilliard equation for a specified number of time steps of size dt.
% 
%INPUTS
    %phi0 = Initial field of chemical states in the domain, created by ch_initialization.
%
%NAME-VALUE PAIRS
    %t_iter = Number of time steps simulated. Default = 1e3.
    %dt = Time step. Default = **2.5e-5 characteristic times** **OR FRACTION OF dt?**.
    %solver_iter = Number of solver iterations per time step. Default = 1e4.
    %tol = Solver tolerance per time step. Default = **1e-6** **OR SCALED BY nx?**.
    %m = Number of mesh points over which the interface exists. Default = 4.
    %boundary = Boundary conditions for the simulation:
    %   'periodic' (default) - flux on one domain border equals negative flux on the opposite border
    %   'neumann' - zero flux on the domain borders
    %c_relax = Number of smoothing relaxations done at the start and end
    %   of each V-cycle. Default = 2 **THIS HAS ALWAYS BEEN TRUE, CORRECT?**;
    %domain = Vector of rightmost and leftmost grid points in x and y.
    %   Format: [xright xleft yright yleft]. Default = [1 0 1 0];
%
%OUTPUT
    %phi_t = Multidimensional array of phi over dt time steps.
    %mass_t = Vector of total mass over dt time steps.
    %E_t = Vector of total energy over dt time steps.

%% Set option defaults and parse inputs

%Set parameter defaults
default_t_iter = 1e3;
default_dt = 2.5e-5;
default_solver_iter = 1e4;
default_tol = 1e-6;
default_m = 4;
default_boundary = 'periodic';
default_c_relax = 2;
default_domain = [1 0 1 0];
default_ns = 10;

CahnHilliard_NMG_parser = inputParser;

%Set general criteria for inputs and name-value pairs
valid_matrix = @(x) ismatrix(x);
valid_integer = @(x) x-floor(x) == 0;
valid_pos_num = @(x) isnumeric(x) && (x > 0);
valid_boundary_type = @(x) strcmpi(x,'periodic') || strcmpi(x,'neumann');
valid_domain_vector = @(x) length(x) == 4;

%Set parser options and valid input criteria
addRequired(CahnHilliard_NMG_parser,'phi0',valid_matrix);
   
addOptional(CahnHilliard_NMG_parser,'t_iter',default_t_iter,valid_integer);
addOptional(CahnHilliard_NMG_parser,'dt',default_dt,valid_pos_num);
addOptional(CahnHilliard_NMG_parser,'solver_iter',default_solver_iter,valid_integer);
addOptional(CahnHilliard_NMG_parser,'tol',default_tol,valid_pos_num);
addOptional(CahnHilliard_NMG_parser,'m',default_m,valid_integer);
addOptional(CahnHilliard_NMG_parser,'boundary',default_boundary,valid_boundary_type);
addOptional(CahnHilliard_NMG_parser,'c_relax',default_c_relax,valid_integer);
addOptional(CahnHilliard_NMG_parser,'domain',default_domain,valid_domain_vector);
addOptional(CahnHilliard_NMG_parser,'write_phi',true);
addOptional(CahnHilliard_NMG_parser,'FileName', "");
addOptional(CahnHilliard_NMG_parser,'ns', default_ns, valid_pos_num);


parse(CahnHilliard_NMG_parser,phi0,varargin{:});

%Extract parsed inputs
phi0 = CahnHilliard_NMG_parser.Results.phi0;
t_iter = CahnHilliard_NMG_parser.Results.t_iter;
dt = CahnHilliard_NMG_parser.Results.dt;
solver_iter = CahnHilliard_NMG_parser.Results.solver_iter;
tol = CahnHilliard_NMG_parser.Results.tol;
m = CahnHilliard_NMG_parser.Results.m;
boundary = CahnHilliard_NMG_parser.Results.boundary;
c_relax = CahnHilliard_NMG_parser.Results.c_relax;
xright = CahnHilliard_NMG_parser.Results.domain(1);
xleft = CahnHilliard_NMG_parser.Results.domain(2);
yright = CahnHilliard_NMG_parser.Results.domain(3);
yleft = CahnHilliard_NMG_parser.Results.domain(4);
write_phi = CahnHilliard_NMG_parser.Results.write_phi;
FileName = CahnHilliard_NMG_parser.Results.FileName;
ns = CahnHilliard_NMG_parser.Results.ns;

%% Define and initialize key simulation parameters

[nx,ny] = size(phi0); %Define number of grid points in x and y
%Define maximum possible V-cycle levels as the largest power of two
%shared by both x and y dimensions
i = 1;
while i < min([floor(log2(nx)) floor(log2(ny))]) && mod(nx,2^i) == 0 ...
    && mod(ny,2^i) == 0
    i = i+1;
end
n_level = i; clear i;
h2 = ((xright-xleft)/nx)*((yright-yleft)/ny); %Define mesh size
epsilon2 = h2*m^2/(2*sqrt(2)*atanh(0.9))^2; %Define ϵ^2
mu = zeros(nx,ny); %Initialize chemical potential
phi_old = phi0; %Initialize prior chemical state
phi_new = phi0; %Initialize next chemical state
% phi_t = zeros(nx,ny,t_iter); %Initialize outputs
mass_t = zeros(t_iter,1);
E_t = zeros(t_iter,1);

%% Run nonlinear multigrid solver
if write_phi
    writematrix(phi_new,sprintf('%s_phi.csv', FileName)); %write IC to file
end
for i = 1:t_iter
    % phi_t(:,:,i) = phi_new;
    mass_t(i) = sum(sum(phi_new))/(h2*nx*ny);
    E_t(i) = discrete_energy(phi_new,h2,nx,ny,epsilon2);
    phi_new = nmg_solver(phi_old,phi_new,mu,nx,ny, ...
        xright,xleft,yright,yleft,c_relax,dt,epsilon2,n_level, ...
        solver_iter,tol,boundary);
    phi_old = phi_new;
    if write_phi
        if rem(i, ns) == 0
            writematrix(phi_new,sprintf('%s_phi.csv', FileName),'WriteMode','append');
        end
    end 
    if mod(i/t_iter*100,5) == 0
        fprintf('%3.0f percent complete\n',i/t_iter*100)
    end
end

final_phi = phi_new;
%Normalize energy to t == 0
E_t = E_t/E_t(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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