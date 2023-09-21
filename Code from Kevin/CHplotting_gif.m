clear
close all hidden
warning off

phi = readmatrix('/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Code from Kevin/phi_spinodal.m','FileType','text');
phidims = size(phi);
phidims(3) = phidims(1)/phidims(2); % Determine number of frames captured
phidims(1) = phidims(2); % Determine size of square grid
phi = reshape(phi, phidims(1), phidims(3), phidims(2)); % Reshape multidimensional array
phi = shiftdim(phi, 2); % Shift dimensions to move frames to the third dimension

% Create a directory to store the individual frames (images)
outputDir = '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C_output/gif_frames';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Create a cell array to store frames
frames = cell(1, phidims(3));

for i = 1:phidims(3)
    phi(:,:,i) = phi(:,:,i); % Transpose so that rows are in the first dimension and columns in the second
    surf(phi(:,:,i),'EdgeColor','none');
    view(2);
    axis([0 size(phi,1) 0 size(phi,2)]);
    colorbar;
    
    % Capture the current frame
    frames{i} = getframe;
    
    % Save the current frame as an image
    frameFilename = fullfile(outputDir, sprintf('frame_%04d.png', i));
    imwrite(frames{i}.cdata, frameFilename);
end

% Define the output GIF filename
gifFilename = '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C_output/phi_spinodal_short_movie.gif';

% Create the GIF from the captured frames
for i = 1:phidims(3)
    im = frames{i}.cdata;
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% Clean up: Delete the individual frame images
delete(fullfile(outputDir, '*.png'));
rmdir(outputDir);

% Close the figure
close;