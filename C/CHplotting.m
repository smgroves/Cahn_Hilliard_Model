clear
close all hidden
warning off

% phi = readmatrix('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Code from Kevin/phi_CPC.m','FileType','text');

phi = readmatrix('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C/phi_CPC_big128.m','FileType','text');
% phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
phidims = size(phi);
phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
phidims(1) = phidims(2); %Determine size of square grid
% phidims(1) = 128;
% phidims(2) = 256;
% phidims(3) = 201;
phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension
myfig = figure();
% hold on
phi_movie = VideoWriter('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/movie_outputs/phi_CPC_big128.mp4','MPEG-4');
phi_movie.Quality = 100; %Highest quality video
open(phi_movie);
% Set consistent color axis limits for the entire movie
clim([-1, 1]);
for i = 1:phidims(3)
    % phi(:,:,i) = phi(:,:,i).'; %Transpose so that rows are in the first dimension and columns in the second
    surf(phi(:,:,i),'EdgeColor','none');
    view(2);
    axis([0 size(phi,1) 0 size(phi,2)]);
    % colormap(Default);
    colorbar;
    % Set the color axis limits for the current frame
    clim([-1, 1]);
    writeVideo(phi_movie,getframe);
end
for i = 1:30
    writeVideo(phi_movie,getframe); %Pause on end frame
end
% surf(double(phi(:,:,size(phi,3))>0),'EdgeColor','none');
% view(2);
% axis([0 size(phi,1) 0 size(phi,2)]);
% colormap(gray);
% colorbar;
% writeVideo(phi_movie,getframe);
% for i = 1:30
%     writeVideo(phi_movie,getframe); %Pause on discretized frame
% end
close(phi_movie);