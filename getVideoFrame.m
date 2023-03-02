function [grayscale,RGB,gradients,options] = getVideoFrame(varargin)
% frame load and processing loop
% this function takes a video input from initializeVideo and grabs a
% specific frame, frame. options for processing are contained in options,
% which is built by setLKoptions;
%
% this function optionally returns the following information:
%
% grayscale : a grayscale image of the frame given by frame
% RGB       : an RGB image of the frame given by frame
% gradients : gradients in X and Y for the frame given by frame
%
% note all outputs are given as scatteredinterpolants and not as gridded data.
% please see matlab doc scatteredInterpolant for more information

video    = varargin{1}; % contains the video structure output by initializeVideo
frame    = varargin{2}; % numeric of the frame of interest to get
options  = varargin{3}; % video processing options
if nargin == 4 % check if user specified to get gradients
    getGradients = varargin{4}; % if they did then use their option
else
    getGradients = false; % otherwise don't take gradients
end

%% switch the different types of formats we have to load the initial frame
% see intializeVideo for more information on formats
switch video.format 
    case 'folder' % folder of images
        RGB = imread(video.File{frame});
    case 'tiff' % tiff stack
        RGB = imread(video.VideoObj,'index',frame);
    case 'movie' % movie file
        RGB = read(video.VideoObj,frame); % current frame
    case {'3dmat' 'mat'}
        RGB = matRead(video,frame);
    case 'bioformats'
        RGB = bfread(video,frame);
    otherwise
        error('unknown input type');
end

%% process the image
if ~options.is3d
    RGB = processRGB(RGB,options); % process the RGB information if desired
    
    grayscaleIM = convertRGBtoGrayscale(RGB,options); % convert the image to grayscale
    % check if user has turned on grayscale image filtering
    if options.ImageProcessing.Filter.on && options.ImageProcessing.Filter.isGray;
        % if it's on process the image according to the given function.
        % filtering can be complex so filtering is specified using a function
        % given in options
        fcn = options.ImageProcessing.Filter.Function; % get the user specified function
        grayscaleIM = feval(fcn,grayscaleIM); % evaluate the given function
    end
else
    switch options.testing.noiseType
        case {'localvar' 'salt & pepper' 'speckle'}
            RGB = imnoise(uint8(RGB),options.testing.noiseType,options.testing.noiseLevel);
        case {'gaussian'}
            RGB = imnoise(uint8(RGB),options.testing.noiseType,0,options.testing.noiseLevel);
        case {'poisson'}
            RGB = imnoise(uint8(RGB),options.testing.noiseType);
    end
    if options.ImageProcessing.FrameResize.Ratio ~= 1 % do we want to resize the images?
        method = options.ImageProcessing.FrameResize.Method; % get the method
        ratio  = options.ImageProcessing.FrameResize.Ratio; % get the ratio
        RGB = ndresize(RGB,ratio,method); % resize the image
    end        
%     if options.testing.addnoise>0 % add random noise to image?
%         sz = size(RGB); % if so get the size of hte random matrix
%         range = (max(RGB(:))-min(RGB(:)));
%         noise = rand(sz)*range*options.testing.addnoise; %generate random noise as function of percent range  
%         RGB = RGB+round(noise); % add the noise to the image
%     end
    if strcmp(video.format,'bioformats')
        grayscaleIM = double(RGB(:,:,:,video.channel));
    else
        grayscaleIM = double(RGB);
    end
end
%% normalize and build the interpolant for the image
if strcmp(options.CameraCalibration,'off') %if we don't have a camera calibration yet
    options.CameraCalibration = calibrateImageCoordinates(grayscaleIM); % make one
end
% then normalize the image
grayscale = normalizeImage(grayscaleIM,options);

%% get the gradients of the image
if getGradients
    %sigma = cellfun(@(x) mean(diff(x)),grayscale.GridVectors);
    sigma = 1./(options.Meshing.boxSize-ones(1,numel(options.Meshing.boxSize)));
    %sigma = [1 1 1];
    if options.is3d
        [Gy,Gx,Gz] = gradient(grayscaleIM,sigma(1),sigma(2),sigma(3));
    else
        [Gy,Gx] = gradient(grayscaleIM,sigma(1),sigma(2)); % gradient returns the coordinates backwards!
    end
    gradients.x = grayscale; % copy the interpolant from above (for speed)
    gradients.y = grayscale; % "" 
    gradients.x.Values = Gx; % replace the values with the gradient
    gradients.y.Values = Gy; % "" 
    if options.is3d
        gradients.z = grayscale; % replace the values with the gradient
        gradients.z.Values = Gz; % copy the interpolant from above (for speed)
    end
else
    gradients = []; % if don't need gradients just output dummy variable
end

end
















%% RGB image processing
function RGB = processRGB(RGB,options)
% check if we need to demosaic the image (some cameras need this)
if options.ImageProcessing.demosaic.on
    % if so demosaic according to the given pattern
    RGB = demosaic(RGB(:,:,1),options.ImageProcessing.demosaic.format);
end
% check if we want to resize the image
if options.ImageProcessing.FrameResize.Ratio ~= 1 % do we want to resize the images?
    method = options.ImageProcessing.FrameResize.Method; % get the method
    ratio  = options.ImageProcessing.FrameResize.Ratio; % get the ratio
    RGB = imresize(RGB,ratio,method); % resize the image
end

% check if user has turned on image filtering
if options.ImageProcessing.Filter.on && ~options.ImageProcessing.Filter.isGray 
    % if it's on process the image according to the given function. 
    % filtering can be complex so filtering is specified using a function 
    % given in options 
    fcn = options.ImageProcessing.Filter.Function; % get the user specified function 
    RGB = feval(fcn,RGB); % evaluate the given function 
end

end

%% conversion to grayscale double for analysis
function grayscale = convertRGBtoGrayscale(RGB,options)
isRGB = ndims(RGB)==3; % check if our input data is actually RGB
if isRGB
    switch options.ImageProcessing.RGB.mode
        case 'builtin' % just convert the image to grayscale using built in command
            grayscale = rgb2gray(RGB);
        case 'mean' % take the mean of the three colors
            grayscale = mean(RGB,3);
        case 'color' % choose a specific color and only take that
            grayscale = RGB(:,:,options.ImageProcessing.RGB.color);
        case 'weighted' % use a weighted average of the RGB values (to do for later)
            weights = options.ImageProcessing.RGB.color;
            grayscale = imapplymatrix(weights, RGB, class(RGB));
    end
    grayscale = double(grayscale); % convert the output to double
else
    grayscale = double(RGB);
end
end

% I don't think I'll be using interpolants anymore. leaving this code now
% incase I come back to it
function [interpolant,K] = normalizeImage(image,options)
% UGH interpolant is slow on massive datasets (images);
% will have to use interp2 to make it reasonably fast.
%
% check if the input is from the computer vision toolbox or a matrix
% if so convert it to a K.fwd, K.inv struct
calibration = options.CameraCalibration;
scatteredData = false;
if isa(calibration,'cameraParameters') || ~isstruct(calibration)
    K = calibrationParameterstoMatrix(calibration); 
elseif isstruct(calibration) % if it's in the right format use it as K
    K = calibration;
end
% check the calibration
if K.fwd(4) ~= 0 && options.ImageProcessing.forceNDGrid
    if options.ImageProcessing.InterpolantWarnings
        warning('Skew enabled in calibration matrix. Forcing off for speed, set options.ImageProcessing.forceNDGrid = false to utilize full scattered data.');
    end
    K.fwd(4) = 0;
elseif ~options.ImageProcessing.forceNDGrid
    if options.ImageProcessing.InterpolantWarnings
        warning('Scattered Data enabled. Please set options.ImageProcessing.forceNDGrid = true for speed.');
    end    
    scatteredData = true;
end
%% need to provide a check in here to force gridded (or a warning!)
sz = size(image); % get the image size
% mesh the grid and linearize the coordinate
if options.is3d
    [X,Y,Z] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
    [x,y,z] = transform3D(X,Y,Z,@(X) K.fwd*X);
    %% create the proper interpolant
    if scatteredData % NOT RECOMMENDED FOR SPEED!!!
        interpolant = scatteredInterpolant(x(:),y(:),z(:),image(:));
    else % gridded interpolants are 2 orders of magnitude faster than scattered data!
        interpolant = griddedInterpolant(x,y,z,image,'cubic');
    end
    
else
    [X,Y] = ndgrid(1:sz(1),1:sz(2));
    [x,y] = transform2D(X,Y,@(X) K.fwd*X);
    %% create the proper interpolant
    
    if scatteredData % NOT RECOMMENDED FOR SPEED!!!
        interpolant = scatteredInterpolant(x(:),y(:),image(:));
    else % gridded interpolants are 2 orders of magnitude faster than scattered data!
        interpolant = griddedInterpolant(x,y,image);
    end
end

end

function RGB = matRead(video,frame)
is3d = numel(video.size)==4;
if isfield(video,'Batches') % we have to process this video in chunks
    global globalVidData %#ok<TLEV> % get or create the global video data variable
    CurrentBatch = video.FrameToBatch(2,video.FrameToBatch(1,:)==frame);
    BatchFrame = video.FrameToBatch(3,video.FrameToBatch(1,:)==frame);
    % check if batch is initialized. if not, find which chunk we need to load first
    if ~globalVidData.BatchInitialized
        % read the proper chunk
        tic
        fprintf('Reading 1st chunk of data ... ');
        batchFrames = video.Batches{CurrentBatch};
        if is3d
            globalVidData.data = video.VideoObj.(video.dataVarName)(:,:,:,batchFrames);
        else
            globalVidData.data = video.VideoObj.(video.dataVarName)(:,:,batchFrames);
        end
        t = toc;
        fprintf('Done! %ds Elapsed\n',round(t))
        globalVidData.Batch = CurrentBatch;
        globalVidData.BatchInitialized = true;
    end
    % check if we're in the right chunk
    if globalVidData.Batch ~= CurrentBatch
        fprintf('Reading %d/%d chunk of data ...... ',CurrentBatch,numel(video.Batches));
        tic
        batchFrames = video.Batches{CurrentBatch};
        if is3d
            globalVidData.data = video.VideoObj.(video.dataVarName)(:,:,:,batchFrames);
        else
            globalVidData.data = video.VideoObj.(video.dataVarName)(:,:,batchFrames);
        end
        globalVidData.Batch = CurrentBatch;
        t = toc;
        fprintf('Done! %ds Elapsed\n',round(t))
    end
    if is3d
        RGB = globalVidData.data(:,:,:,BatchFrame);
    else
        RGB = globalVidData.data(:,:,BatchFrame);
    end
else
    if is3d
        RGB = video.VideoObj.(video.dataVarName)(:,:,:,frame);
    else
        RGB = video.VideoObj.(video.dataVarName)(:,:,frame);
    end
end
if ~isa(RGB,'double') 
    RGB = double(RGB);
end
end

function K = calibrateImageCoordinates(image)
    nd = ndims(image);
    sz = size(image);
    if nd == 3
        xv = 1:sz(1);
        yv = 1:sz(2);
        zv = 1:sz(3);
        [X.x,X.y,X.z] = ndgrid(xv,yv,zv);
    elseif nd ==2 
        xv = 1:sz(1);
        yv = 1:sz(2);
        [X.x,X.y] = ndgrid(xv,yv);
    end
    K = normalizeCoordinateSystem(X);
    %K.fwd = Kt.inv;
    %K.inv = Kt.fwd;
end

function RGB = bfread(video,frame)
RGB = zeros([video.size numel(video.channels)],'like',(video.channels{1}));
for i = 1:numel(video.channels)
    RGB(:,:,:,i) = video.channels{i}(:,:,:,frame);
end

if ~isa(RGB,'double')
    RGB = double(RGB);
end
end

function Aresize = ndresize(A,ratio,method)
if strcmp(method,'lanczos3') % method was left default, use cubic
    method = 'cubic';
end
sz = size(A);
newSize = round(ratio*sz);
xv = linspace(1,sz(1),newSize(1));
yv = linspace(1,sz(2),newSize(2));
zv = linspace(1,sz(3),newSize(3));
[x,y,z] = ndgrid(xv,yv,zv);
Ainterp = griddedInterpolant(A,method);
Aresize = Ainterp(x,y,z);
end