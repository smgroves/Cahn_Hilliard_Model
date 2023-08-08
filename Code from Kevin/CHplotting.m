clear
close all hidden
warning off

phi = readmatrix('/Users/kaj5f/Desktop/Cahn-Hilliard solver Mathematics 2020/phi.m','FileType','text');
phidims = size(phi);
phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
phidims(1) = phidims(2); %Determine size of square grid
phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension
phi_movie = VideoWriter('/Users/kaj5f/Desktop/phi_movie.mp4','MPEG-4');
phi_movie.Quality = 100; %Highest quality video
open(phi_movie);
for i = 1:phidims(3)
    phi(:,:,i) = phi(:,:,i); %Transpose so that rows are in the first dimension and columns in the second
    surf(phi(:,:,i),'EdgeColor','none');
    view(2);
    axis([0 size(phi,1) 0 size(phi,2)]);
    colormap(viridis);
    colorbar;
    writeVideo(phi_movie,getframe);
end
for i = 1:30
    writeVideo(phi_movie,getframe); %Pause on end frame
end
surf(double(phi(:,:,size(phi,3))>0),'EdgeColor','none');
view(2);
axis([0 size(phi,1) 0 size(phi,2)]);
colormap(gray);
colorbar;
writeVideo(phi_movie,getframe);
for i = 1:30
    writeVideo(phi_movie,getframe); %Pause on discretized frame
end
close(phi_movie);