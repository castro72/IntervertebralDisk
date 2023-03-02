function options = setLKoptions(varargin)
%% General Options
% Set up the warp type we're using. Please specify functions in warpingFunctions!
options.WarpType = 'affine3d'; % currently fully implemented: 'affine2d' 'affine3d' 'affine1d' 'affine1dmanifold'
%% Image Processing
options.ImageProcessing.demosaic.on = false; % is the image in bayer format
options.ImageProcessing.demosaic.format = 'gbrg'; % bayer method decode
options.ImageProcessing.FrameResize.Method = 'lanczos3'; % averaging method to resize the frame
options.ImageProcessing.FrameResize.Ratio = 1 ; % frame resize ratio
options.ImageProcessing.Filter.on = false ; % apply an additional image filter
options.ImageProcessing.Filter.Function = ''; % filter function
options.ImageProcessing.Filter.isGray = 0; % is the filter for grayscale only?
options.ImageProcessing.RGB.mode = 'builtin'; % RGB conversion mode ('builtin' 'mean' 'color' 'weighted')
options.ImageProcessing.RGB.color = []; % color to specify for choosing an RGB color, or specify as a 3 vector of weights
options.ImageProcessing.forceNDGrid = true; % highly recommended for speed. turn off if you want to allow non-rectangular coordinates
options.ImageProcessing.InterpolantWarnings = true; % highly recommended for speed. turn off if you want to allow non-rectangular coordinates

%% Testing
options.testing.noiseLevel = 0; % add random noise to the image?
options.testing.noiseType = 'none'; % {'gaussian'  'poisson' 'salt & pepper' 'speckle'}
options.nxc = false; % do we want to run a normal cross correlation study in parrallel with the LK solver?
%% Meshing
options.Meshing.ConsistentDimensions = true; % set true if all analysis regions are the same size
options.Meshing.mode = 'boxes';
options.Meshing.spacing = [40 40 40]; % how far apart should boxes be spaced
options.Meshing.boxSize = [51 51 51]; % should be odd
options.Meshing.AnalysisRegionPrecentage = 0.50; % relative size of a box about the center to automatically mesh
options.Meshing.Mask = []; % mask these points out
%% Optimization
options.Optimization.maximumIterations = 101; % max number of iterations to perform
options.Optimization.epsilon = 0.001; % epsilon that must be reached before optimzation is stopped
options.Optimization.LM = false; % Enable or disable Levenberg-Marquardt optimization.
options.Optimization.recovery = false; % enable region recovery functions
options.Optimization.minMatrixCondition = 0.5; % minimum matrix condition to allow before deciding the optimizaton has failed.
%% Mechanics
options.mechanics.LSF_Points = 5;
options.mechanics.fixDim = true;
%% Plotting
options.plotting.HeatMapAlpha = .2; 

%% General Options
options.statusBars = true;
%% Processing based on other options
% check if we're using a meshing method with consistent dimensions
if strcmp(options.Meshing.mode,'boxes')        
    options.Meshing.ConsistentDimensions = true;
end

%% processing name/value pairs as inputs, code based on: 
% http://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('need propertyName/propertyValue pairs')
end

% Camera Calibration specific arguments
inNames = varargin(1:2:nArgs); % get the input names
% see if CameraCalibration is one of them
calibLoc = strcmp(inNames,'CameraCalibration'); 
if any(calibLoc) %checks to see if it is true (1)
    calibLoc = find(calibLoc);
    calibrationSession = varargin{2*calibLoc};
    varargin([2*calibLoc 2*calibLoc-1]) = [];
else
    isValid = false;
    while ~isValid
        [calibFile,calibPath] = uigetfile('*.mat','Please select a camera calibration','Camera Calibration File');
        load([calibPath calibFile ]); % load the session
        % check if valid
        isValid = exist('calibrationSession','var');
        if ~isValid
            warning([calibFil ' is not a calibration session, please reselect a calibration session']);
        end
    end
end
if strcmp(calibrationSession,'off')
    options.CameraCalibration = 'off';
else
    options.CameraCalibration = calibrationSession.CameraParameters;
end
% read the acceptable names
optionNames = getAllFieldNames(options);

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = pair{1}; % make case insensitive

   if any(strcmp(inpName,optionNames))
      % overwrite options. If you want you can test for the right class here
      % Also, if you find out that there is an option you keep getting wrong,
      % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options = replacePair(inpName,pair{2},options);
      %options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

% check if we are using 3d
options.is3d = check3d(options);
% get the warping functions
options = warpingFunctions(options);
end

function [allFields,newFields]= getAllFieldNames(struct,parent,allFields)
% recursive function to get all fields and their children of a structure.

if nargin == 1; % if we're on the first try
    parent = ''; % there is no parent
    allFields = {}; % and no fields yet
else
    parent = [parent '.']; % otherwise add the dot to the parent
end

if isstruct(struct) % it's a structure, process it and get the names
    CurrentFieldNames = fieldnames(struct); % get the field names
    for i = 1:numel(CurrentFieldNames) % loop through each one
        name = CurrentFieldNames{i}; % get the name
        [allFields,newFields] = getAllFieldNames(struct.(name),[parent name],allFields); % process it
        if numel(newFields) ~= 0  % if there are new fields
            newFields = cellfun(@(x) [parent '.' x],newFields,'UniformOutput',false); % add the parent to the field
        end
    end
else  % if it's not a structure, add the field to the result
    allFields = [allFields parent(1:end-1)]; % add the result (removing the period added above)
    newFields = []; % no new fields
end
end

function options = replacePair(inpName,pair,options)
% replace the value in the options with the new value
dots = strfind(inpName,'.'); % get the structure
dots = [1 dots numel(inpName)]; % find the dots
structs = {}; % make a dummy variable to loop over
for d = 2:numel(dots)
    structs = [structs inpName(dots(d-1):dots(d))]; %#ok<AGROW> % add each value
end
structs = cellfun(@(x) regexprep(x,'\.',''),structs,'UniformOutput',false); % get rid of the dots
options = setfield(options,structs{:},pair); % set the field value
end

function is3d = check3d(options)
% function to check if the LK method being used is 3d
warpType = options.WarpType; % get our warp type
valid3dWarps = {'affine3d' 'projective3d' 'affine3dPrinciple'}; % check the valid3d warp types
is3d = any(strcmp(warpType,valid3dWarps));
end