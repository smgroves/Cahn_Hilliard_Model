function ch_movie(t_out,phi_t,varargin)
%This function creates a red-white-blue video trajectory of
%chemical states in the current directory
%
%INPUTS
% t_out = Vector of output times corresponding to each frame in phi_t.
% phi_t = Multidimensional array of chemical states at the corresponding times in t_out.
%
%NAME-VALUE PAIRS
% filename: String for the output movie filename. Default = 'ch_movie'.
% filetype: Movie file type compatible with VideoWriter. Defaults to 'MPEG-4'.
%
%OUTPUT
% None except for the saved red-white-blue video.

warning off

%% Set option defaults and parse inputs

%Set parameter defaults
default_filename = 'ch_movie';
default_filetype = 'MPEG-4';

ch_movie_parser = inputParser;

%Set general criteria for inputs and name-value pairs
valid_3darray = @(x) ismatrix(x(:,:,1));
valid_filename = @(x) ischar(x) || isstring(x);
valid_filetype = @(x) strcmpi(x,'MPEG-4') || strcmpi(x,'Motion JPEG AVI') ...
    || strcmpi(x,'Uncompressed AVI');
valid_time = @(x) isvector(x) && ~isempty(x);

%Set parser options and valid input criteria
addRequired(ch_movie_parser,'t_out',valid_time);
addRequired(ch_movie_parser,'phi_t',valid_3darray);
addParameter(ch_movie_parser,'filename',default_filename,valid_filename);
addParameter(ch_movie_parser,'filetype',default_filetype,valid_filetype);

parse(ch_movie_parser,t_out,phi_t,varargin{:});

%Extract parsed inputs
t_out = ch_movie_parser.Results.t_out;
phi_t = ch_movie_parser.Results.phi_t;
filename = ch_movie_parser.Results.filename;
filetype = ch_movie_parser.Results.filetype;

% Check that length of t_out matches the number of frames
if length(t_out) ~= size(phi_t,3)
    error('Length of t_out must match the number of frames in phi_t (size(phi_t,3)).');
end

phi_movie = VideoWriter(strcat(cd,'/',filename),filetype); %Extension will be automatically appended
if strcmpi(filetype,'MPEG-4') || strcmpi(filetype,'Motion JPEG AVI')
    phi_movie.Quality = 100; %Highest quality video compression
end
open(phi_movie);

% Choose a format specification
formatSpec = '%10.6f'; 

for i = 1:size(phi_t,3)
    image(transpose(phi_t(:,:,i)),'CDataMapping','scaled');
    colorbar; 
    axis square;
    clim([min(phi_t(:,:,1),[],'all') max(phi_t(:,:,1),[],'all')]); %Set color axis to initial conditions
    g = gca;
    %Scale and center display roughly to the size of the mesh,
    %with extra space horizontally for the color bar
    % g.Position = [0.5-size(phi_t,2)/1000 0.5-size(phi_t,1)/1000 ...
    %     2.2*size(phi_t,2)/1000 2*size(phi_t,1)/1000];
    axis([0 size(phi_t,2) 0 size(phi_t,1)]);
    colormap(redbluecmap);
    colorbar;
    
    % Add a title showing the current time
    % Use sprintf to ensure uniform length and alignment
    timeStr = sprintf(formatSpec, t_out(i));
    title(['Time = ' timeStr]);

    frame = getframe(gcf);
    writeVideo(phi_movie,frame);
end

close(phi_movie);
