% Assume 'data' is a 3D matrix where each slice data(:,:,t) is the data at time point t
% Replace with actual loading mechanism
% data = load('your_data_file.mat'); % Load your actual data here

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry";
total_time=.025;
dt_name="2.5e-5";
dt = str2double(dt_name);
ns = 10;
dt_in_movie = dt*ns;
timesteps=total_time/dt;
name=sprintf('phi_256_%s_1.0e-5__CPC_10_cohesin_4_eps_0.007504684956431058',string(timesteps))
dtout=10;
phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
phidims = size(phi);
phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
phidims(1) = phidims(2); %Determine size of square grid
phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension
data = phi;

colors = {
    "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
    "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
};
% Assume 'data' is a 3D matrix where each slice data(:,:,t) is the data at time point t
% Replace with actual loading mechanism
% data = load('your_data_file.mat'); % Load your actual data here

% Initialize parameters
numTimePoints = size(data, 3); % Number of time points
dropletData = struct(); % Structure to store droplet data
dt = 0.1; % Time interval (for example, 0.1 seconds)

% Loop over each time point
for t = 1:numTimePoints
    % Extract the data for the current time point
    currentData = data(:,:,t);
    
    % Find the 0-level contour
    contourMatrix = contourc(currentData, [0, 0]);
    
    % Initialize a list to store droplets for the current time point
    droplets = [];
    startIdx = 1;

    % Parse the contour matrix to extract droplet contours
    while startIdx < size(contourMatrix, 2)
        % Read the number of points in the current contour line
        numPoints = contourMatrix(2, startIdx);
        dropletContour = contourMatrix(:, startIdx+1:startIdx+numPoints);
        
        % Calculate the center and radius of the droplet
        centerX = mean(dropletContour(1, :));
        centerY = mean(dropletContour(2, :));
        distances = sqrt((dropletContour(1, :) - centerX).^2 + (dropletContour(2, :) - centerY).^2);
        radius = mean(distances);
        
        % Store droplet information
        droplet = struct('center', [centerX, centerY], 'radius', radius);
        droplets = [droplets, droplet]; %#ok<AGROW>
        
        % Move to the next contour
        startIdx = startIdx + numPoints + 1;
    end
    
    % Store the droplet data for the current time point
    dropletData(t).droplets = droplets;
end

% Initialize a structure array to track droplets across time points
trackedDroplets = struct('radius', {}, 'time', {}, 'center', {}, 'id', {});

% Define an appropriate threshold for droplet matching
some_threshold = 10; % This value should be adjusted based on your data

% Track the droplets
for t = 1:numTimePoints
    currentDroplets = dropletData(t).droplets;
    
    if t == 1
        % Initialize tracked droplets with the first time point's data
        for i = 1:length(currentDroplets)
            trackedDroplets(i).radius = currentDroplets(i).radius;
            trackedDroplets(i).time = t * dt; % Store the actual time
            trackedDroplets(i).center = currentDroplets(i).center;
            trackedDroplets(i).id = i; % Assign a unique ID
        end
    else
        % Extract previous droplet centers
        previousCenters = arrayfun(@(d) d.center, trackedDroplets, 'UniformOutput', false);
        previousCenters = cell2mat(previousCenters(:)'); % Convert to matrix

        % Extract current droplet centers
        currentCenters = arrayfun(@(d) d.center, currentDroplets, 'UniformOutput', false);
        currentCenters = cell2mat(currentCenters(:)'); % Convert to matrix

        % Calculate the distance matrix between previous and current centers
        distMatrix = pdist2(previousCenters, currentCenters);

        % Assign current droplets to the nearest previous droplet
        [minDists, minIdxs] = min(distMatrix, [], 1);
        
        matched = false(1, length(currentDroplets));
        
        for j = 1:length(currentDroplets)
            if minDists(j) < some_threshold
                % Update the tracked droplet
                idx = minIdxs(j);
                trackedDroplets(idx).radius = [trackedDroplets(idx).radius, currentDroplets(j).radius];
                trackedDroplets(idx).time = [trackedDroplets(idx).time, t * dt]; % Store the actual time
                trackedDroplets(idx).center = [trackedDroplets(idx).center; currentDroplets(j).center];
                matched(j) = true;
            end
        end
        
        % Add unmatched droplets as new droplets
        for j = 1:length(currentDroplets)
            if ~matched(j)
                newIdx = length(trackedDroplets) + 1;
                trackedDroplets(newIdx).radius = currentDroplets(j).radius;
                trackedDroplets(newIdx).time = t * dt; % Store the actual time
                trackedDroplets(newIdx).center = currentDroplets(j).center;
                trackedDroplets(newIdx).id = newIdx; % Assign a unique ID
            end
        end
    end
end

% Plot the radius of each droplet over time with custom x-axis labels
figure;
hold on;
for i = 1:length(trackedDroplets)
    plot(trackedDroplets(i).time, trackedDroplets(i).radius, '-o');
end
xlabel('Time (seconds)');
ylabel('Radius');
title('Radius of Each Droplet Over Time');
hold off;
