function phi_new = nmg_solver(phi_old,phi_new,mu,nx,ny, ...
    xright,xleft,yright,yleft,c_relax,dt,epsilon2,n_level,max_iter,tol, ...
    boundary)
%This function uses the nonlinear multigrid method to solve the 
%Cahn-Hilliard equation for the next time step.
%
%INPUTS
    %phi_old = Prior chemical state.
    %phi_new = Next chemical state.
    %mu = Chemical potential
    %nx = Number of grid points in x.
    %ny = Number of grid points in y.
    %xright = Rightmost grid point in x.
    %xleft = Leftmost grid point in x.
    %yright = Rightmost gridpoint in y.
    %yleft = Leftmost gridpoint in y.
    %ilevel = Current V-cycle level
    %c_relax = Number of smoothing relaxations.
    %dt = Time step.
    %epsilon2 = Squared interfacial transition distance (see ϵ^2 term in Lee et al., Mathematics 2020)
    %n_level = Total number of V-cycles.
    %max_iter = Maximum number of iterations for numerical convergence
    %tol = Tolerance for numerical convergence
    %boundary = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %phi_new = Next chemical state.

iter = 0;
resid2 = 1;

%Calculate source terms
s_mu = zeros(nx,ny);
phi_lap = nmg_laplace(phi_old,nx,ny,xright,xleft,yright,yleft,boundary);
s_phi = phi_old/dt - phi_lap;

while resid2 > tol && iter < max_iter
    [phi_new,mu] = nmg_vcycle(phi_new,mu,s_phi,s_mu,nx,ny, ...
        xright,xleft,yright,yleft,1,c_relax,dt,epsilon2,n_level,boundary); %V-cycle starting a level 1
    resid2 = error2(phi_old,phi_new,mu,nx,ny,xright,xleft,yright,yleft,dt,boundary);
    iter = iter+1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res2 = error2(phi_oldf,phi_newf,muf,nxf,nyf,xr,xl,yr,yl,dtf,bc)
%This function computes the 2D residual for phi.
%
%INPUTS
    %phi_oldf = Prior chemical state.
    %phi_newf = Next chemical state.
    %muf = Chemical potential.
    %nxf = Number of grid points in x.
    %nyf = Number of grid points in y.
    %xr = Rightmost grid point in x.
    %xl = Leftmost grid point in x.
    %yr = Rightmost gridpoint in y.
    %yl = Leftmost gridpoint in y.
    %dtf = Time step.
    %bc = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %res2 = Normalized residual

rr = muf - phi_oldf; %Calculate the starting residual
sor = nmg_laplace(rr,nxf,nyf,xr,xl,yr,yl,bc); %Calculate the source term from rr
rr = sor - (phi_newf-phi_oldf)/dtf; %Update the residual
res2 = sqrt(sum(sum(rr.^2))/(nxf*nyf)); %Take Frobenius norm of sum of squared error

end