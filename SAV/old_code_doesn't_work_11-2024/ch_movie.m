function ch_movie(phi_t,varargin)
%This function creates a red-white-blue video trajectory of
%chemical states in the current directory
%
%INPUTS
    %phi_t = Multidimensional array of chemical states.
%
%NAME-VALUE PAIRS
    %filename = Number of movie file to be saved. Default = 'ch_movie'.
    %filetype = Any permissible movie profile in VideoWriter. Most useful:
    %   'MPEG-4' (default) - MPEG-4 file with H.264 encoding; best
    %       compression but may cause artifacts with image fields
    %       of high spatial frequency
    %   'Motion JPEG AVI' - AVI file using Motion JPEG encoding; max
    %       quality JPEG compression
    %   'Uncompressed AVI' - Uncompressed AVI file with RGB24 video;
    %       largest file size
%
%OUTPUT
    %None except for the saved red-white-blue video

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

%Set parser options and valid input criteria
addRequired(ch_movie_parser,'phi_t',valid_3darray);
   
addParameter(ch_movie_parser,'filename',default_filename,valid_filename);
addParameter(ch_movie_parser,'filetype',default_filetype,valid_filetype);

parse(ch_movie_parser,phi_t,varargin{:});

%Extract parsed inputs
phi_t = ch_movie_parser.Results.phi_t;
filename = ch_movie_parser.Results.filename;
filetype = ch_movie_parser.Results.filetype;

phi_movie = VideoWriter(strcat(cd,'/',filename),filetype); %Extension will be automatically appended
if strcmpi(filetype,'MPEG-4') || strcmpi(filetype,'Motion JPEG AVI')
    phi_movie.Quality = 100; %Highest quality video compression
end
open(phi_movie);
for i = 1:size(phi_t,3)
    image(transpose(phi_t(:,:,i)),'CDataMapping','scaled');
    colorbar; 
    axis square;
    clim([min(phi_t(:,:,1),[],'all') max(phi_t(:,:,1),[],'all')]); %Set color axis to initial conditions
    g = gca;
    %Scale and center display roughly to the size of the mesh,
    %with extra space horizontally for the color bar
    g.Position = [0.5-size(phi_t,2)/1000 0.5-size(phi_t,1)/1000 ...
        2.2*size(phi_t,2)/1000 2*size(phi_t,1)/1000];
    axis([0 size(phi_t,2) 0 size(phi_t,1)]);
    colormap(redbluecmap);
    colorbar;
    writeVideo(phi_movie,getframe);
end
close(phi_movie);