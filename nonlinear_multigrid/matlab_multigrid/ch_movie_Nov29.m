function ch_movie(phi_t)
%This function creates a pseudocolore mp4 video of trajectory of
%chemical states
%
%INPUTS
    %phi_t = Multidimensional array of chemical states.
%
%NAME-VALUE PAIRS
    %t_iter = Number of time steps simulated. Default = 1e3.
    %dt = Time step. Default = 2.5e-5 characteristic times.
    %solver_iter = Number of solver iterations per time step. Default = 1e4.
    %tol = Solver tolerance per time step. Default = 1e-6.
    %m = Number of mesh points over which the interface exists. Default = 8.
    %c_relax = Number of smoothing relaxations done at the start and end
    %   of each V-cycle. Default = 2;
    %domain = Vector of rightmost and leftmost grid points in x and y.
    %   Format: [xright xleft yright yleft]. Default = [1 0 1 0];
%
%OUTPUT
    %phi_t = Multidimensional array of phi over dt time steps.
    %mass_t = Vector of total mass over dt time steps.
    %E_t = Vector of total energy over dt time steps.

warning off

phi_movie = VideoWriter(strcat(cd,'/ch_movie.mp4'),'MPEG-4');
phi_movie.Quality = 100; %Highest quality video
open(phi_movie);
for i = 1:size(phi_t,3)
    surf(phi_t(:,:,i),'EdgeColor','none');
    view(2);
    axis([0 size(phi_t,2) 0 size(phi_t,1)]);
    colormap(viridis);
    colorbar;
    writeVideo(phi_movie,getframe);
end
close(phi_movie);