function video = initializeVideo(inputType,getFrameTimes,file)
%
% create movie structure which describes the structure of the video file.
% video.File = file or files pointing to the video file
% video.FrameTimes = time of each frame (in seconds)
% video.format = the format of video
% video.videoObj = video object for matlab (if needed)
% use this function in conjuction with getFrame to load an individual frame

%%
if nargin == 3 % if the file is provided
    switch inputType
        case 'folder'
            video = folderLoad(file,getFrameTimes);
        case 'tiff'
            video = tiffLoad(file,getFrameTimes);
        case 'video'
            video = movieLoad(file,getFrameTimes);
        case {'mat' '3dmat'}
            video = matLoad(file,getFrameTimes);
        case 'bioformats'
            video = bfload(file,getFrameTimes);
    end
else % if nargin == 1, we need to load the file too
    switch inputType % switch the type of input
        case 'folder'
            vidDir = uigetdir;
            video = folderLoad(vidDir,getFrameTimes);
        case 'tiff'
            [vidFile,vidPath] = uigetfile('*.*');
            video = tiffLoad([vidPath vidFile],getFrameTimes);
        case 'video'
            [vidFile,vidPath] = uigetfile('*.*');
            video = movieLoad([vidPath vidFile],getFrameTimes);
    end
end
video.Frames = 1:video.NumberOfFrames;
end

%% below are the loading and processing functions for each video format
function video = folderLoad(movDir,getFrameTimes)
video.format = 'folder'; % set the format
fs = dir(movDir); % get all files
fs([fs.isdir]) = []; % get rid of anything that's not a directory
n = 0; % start counter
for i = 1:numel(fs) % loop through each file
    [~,~,ext] = fileparts(fs(i).name);
    isImage = any(strcmp(ext,{'.jpg','.bmp','.tif'}));
    if isImage
        n = n+1;
        if isunix
            slash = '/';
        elseif ispc
            slash = '\';
        end
        video.File{n} = [movDir slash fs(i).name]; % add it to the list
    end
end
video.NumberOfFrames = numel(video.File); % save the number of frames
if getFrameTimes
    video.FrameTimes = getBMPFrameTimes(video.File); % save the frame times
else 
    video.FrameTimes = 0:video.NumberOfFrames-1;
end
video.VideoObj = 'files'; % placeholder 

end

function video = tiffLoad(vidFile,getFrameTimes)
    % basic function for loading frame from a TIFF
    video.format = 'tiff';
    video.File = vidFile;
    video.NumberOfFrames = numel(imfinfo(vidFile)); % crude but it works
if getFrameTimes
    video.FrameTimes = 0:video.NumberOfFrames-1;
else 
    video.FrameTimes = 0:video.NumberOfFrames-1;
end
    video.FrameTimes = 0:video.NumberOfFrames-1;
    video.VideoObj = 'tiff';
end

function video = movieLoad(vidFile,getFrameTimes)
video.format = 'movie'; % single input movie file
video.File = vidFile; % save the file
video.VideoObj = VideoReader(vidFile); % create the video object
video.NumberOfFrames = video.VideoObj.NumberOfFrames; % get the number of frames
% get FrameTimes (Based off of frame rate and number of frames)
if getFrameTimes
video.FrameTimes = 0 : 1/video.VideoObj.FrameRate : ...
                   video.NumberOfFrames/video.VideoObj.FrameRate; 
else
     video.FrameTimes =  0:video.NumberOfFrames-1;
end
end

function video = matLoad(vidFile,getFrameTimes)
% function written to grab 3d data from a mat file (specific!)
forceDisk = false; % change this to force the program to load from disk
memoryPerecent = .5; % 10% max available memory should be loaded.
video.format = 'mat';
if ischar(vidFile)
    video.File = vidFile;
    video.VideoObj = matfile(vidFile);
    % check the video size and see if we want to load it in memory or load the
    % entire thing
    vidInfo = whos(video.VideoObj);
    vidMemSize = sum([vidInfo.bytes]); % video data size in bytes
    dataidx = find(strcmp('data',{vidInfo.name}));
    if isempty(dataidx)
        if numel({vidInfo.name})>1;
            dataidx = listdlg('PromptString','Please Select Data Variable',...
                'SelectionMode','single','ListString',{vidInfo.name});
            
        else iscell(vidInfo.name)
            dataidx = 1;
        end
    end
    dataVarName = vidInfo(dataidx).name;
%     if isunix
%         [r,w] = unix('free | grep Mem');
%         stats = str2double(regexp(w, '[0-9]*', 'match'));
%         mem.PhysicalMemory.Total = stats(1);
%         mem.PhysicalMemory.Available = (stats(3)+stats(end))*1024;
%     else
%         [~,mem] = memory; % check system memory
%     end
%     if (vidMemSize/mem.PhysicalMemory.Available > memoryPerecent) || forceDisk % the array is more than 10% the available memory, just load it in
%         global globalVidData
%         statsString = sprintf('(%d / %d MB req/free)',round(vidMemSize/1024^2),round(mem.PhysicalMemory.Available/1024^2));
%         warning(['Array is more than half available memory, will load frames in chunks from disk ' statsString]);
%         %% find out the size of a single frame
%         % find the data entry
%         
%         vidSz = vidInfo(dataidx).size;
%         nFrames = vidSz(4);
%         vidClass = vidInfo(dataidx).class;
%         classTypes = {'single' 'double' 'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64'};
%         classBytes = [ 4        8        1      2       3       4       1       2        3        4];
%         vidBytes = classBytes(strcmp(vidClass,classTypes));
%         vidFrameSize = vidBytes * vidSz(1) * vidSz(2) * vidSz(3); % get video frame size in bytes
%         % find out how many frames will take up half of available memory
%         maxFramesLoaded = floor(mem.PhysicalMemory.Available/vidFrameSize*memoryPerecent);
%         nBatches = ceil(nFrames/maxFramesLoaded); % find out how many bathces we should do
%         % break up frames into even batches to load
%         bStart = 1;
%         batches = cell(1,nBatches);
%         FrameToBatch = [1:nFrames ; zeros(1,nFrames); zeros(1,nFrames)];
%         for i = 1:nBatches
%             bEnd = bStart+ceil(nFrames/(nBatches))-1;
%             if bEnd>nFrames
%                 bEnd = nFrames;
%             end
%             batches{i} = bStart:bEnd;
%             FrameToBatch(2,batches{i}) = i;
%             FrameToBatch(3,batches{i}) = 1:numel(batches{i});
%             bStart = bEnd+1;
%         end
%         video.FrameToBatch = FrameToBatch; % give pointers to each
%         video.Batches = batches; % add the batches to the video
%         globalVidData.BatchInitialized = false; % write that the batch isn't initialized to the file
%         video.VideoObj.Properties.Writable = true;
%         video.VideoObj.FrameRate = 0.01; % add a framerate
%         video.VideoObj.Properties.Writable = false;
%         video.size = size(video.VideoObj,dataVarName);
%         video.dataVarName = dataVarName;
%     else
        video.dataVarName = vidInfo(dataidx).name;
        fprintf('Reading entire video into memory...')
        video.VideoObj = load(vidFile);
        video.size = size(video.VideoObj.(video.dataVarName));
        fprintf(' Done!\n')
%     end
else
    video.dataVarName = 'data';
    video.VideoObj.data = vidFile;
    video.size = size(video.VideoObj.(video.dataVarName));
end
% keyboard;
% video.NumberOfFrames = video.size(10);
video.NumberOfFrames = video.size(end);
if numel(video.size) == 3
    video.dims = {'x' 'y' 't'};
else
    video.dims = {'x' 'y' 'z' 't'};
end
if getFrameTimes
    video.FrameTimes = 0 : 1/video.VideoObj.FrameRate : ...
                   video.NumberOfFrames/video.VideoObj.FrameRate; 
else
     video.FrameTimes =  0:video.NumberOfFrames-1;
end
end

function times = getBMPFrameTimes(frames)
% this function checks to see if the input is multiple cells (multiple
% files) and processes all of them. If it's just a single input, process
% the single file name and return the time.
if iscell(frames) % if we have a cell array perform the operation on each cell
    times = cellfun(@(x) getBMPFrameTime(x), frames); % get each frame time for each cell
    times = times-times(1); % subtract the first time to get the relative times
else
    times = getframetime(frames);
end

end

function time = getBMPFrameTime(frame_name)
% function for extracting a single frame time
% getBMPFrameTime pulls the time out of the filename in
% YYYY-MM-DD-HH-MM-SS-MLS, then converts to datenum format and unix
% format (i.e., seconds since Jan 01, 1970).
hyphens = regexp(frame_name,'-');
YYYY    = str2double(frame_name(hyphens(end-5)-4 : hyphens(end-5)-1));
MM      = str2double(frame_name(hyphens(end-4)-2 : hyphens(end-4)-1));
DD      = str2double(frame_name(hyphens(end-3)-2 : hyphens(end-3)-1));
HH      = str2double(frame_name(hyphens(end-2)-2 : hyphens(end-2)-1));
MN      = str2double(frame_name(hyphens(end-1)-2 : hyphens(end-1)-1));
SS      = str2double(frame_name(hyphens(end)  -2 : hyphens(end)  -1));
MLS     = str2double(frame_name(hyphens(end)  +1 : hyphens(end)  +3));

frame_date_M    = datenum(YYYY, MM, DD, HH, MN, SS + 0.001*MLS);
    % datenum outputs days since Jan 0, 0000 (beginning of AD)
time = 86400 * (frame_date_M - datenum('01-Jan-1970'));
    % converts to Unix time = seconds since Jan 1, 1970
end

function video = bfload(file,getFrameTimes) % load a bioformat file

data = bfopen(file); % read the file
strings = data{1,1}(:,2); 
s=textscan(strings{1},'%s plane%s Z=%s C=%s T=%s','Delimiter',';');
ind = 0;
if isempty(s{3})
    s=textscan(strings{1},'%s %s plane%s Z?=%s C?=%s T?=%s','Delimiter',';');
    ind = 1;
end
Zdepth = 1/str2num(s{3+ind}{1});
Ttotal = 1/str2num(s{5+ind}{1});
Ctotal = 1/str2num(s{4+ind}{1});
for i = 1:numel(strings)   
    s=textscan(strings{i},'%s plane%s Z=%s C=%s T=%s','Delimiter',';');
    if isempty(s{3})
        s=textscan(strings{i},'%s %s plane%s Z?=%s C?=%s T?=%s','Delimiter',';');
    end
    Z = str2num(s{3+ind}{1})*Zdepth;
    T = str2num(s{5+ind}{1})*Ttotal;
    C = str2num(s{4+ind}{1})*Ctotal;
    ZTC(:,i) = [Z T C];
end

% check size of image
% preallocate
imsize = size(data{1}{1});
imclass = class(data{1}{1});
channels = cell(1,Ctotal);
cmax = zeros(1,Ctotal);
for Ci = 1:Ctotal
    channels{Ci} = zeros(imsize(1),imsize(2),Zdepth,Ttotal,imclass);
    for Ti = 1:Ttotal
        for Zi = 1:Zdepth
            % find the proper image
            imnum = ismember(round(ZTC'),[Zi Ti Ci],'rows');
            % insert it
            channels{Ci}(:,:,Zi,Ti) = data{1,1}{imnum,1};
        end
    end
    cmax(Ci) = max(channels{Ci}(:));
end

%% gathering the remaining info
bitDepth = log2(double(max(cmax+1)));
bf.channels = channels;
omeMeta = data{1, 4};
bf.voxelSizeX = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in µm
bf.voxelSizeY = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in µm
bf.voxelSizeZ = double(omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM)); % in µm
bf.voxelUnit = 'µm';
bf.bitDepth = bitDepth; 
metadata = data{1, 2};
bf.details = metadata.get('Global sSpecSettings');
bf.NumberOfFrames = Ttotal;


% choose channel
f=figure;
title('Please select a channel to analyze')
colormap gray
for i = 1:Ctotal
    subplot(1,Ctotal,i)
    im(i) = imagesc(bf.channels{i}(:,:,round(Zdepth/2),1));
    axis square
    fcn = @axeschosen;
    set(im(i),'ButtonDownFcn',{fcn,i});
end
%set(f,'OutputFcn','varargout{1} = x;')
uiwait
SelectedChannel = get(f,'UserData');
close(f)
bf.channel = SelectedChannel;
bf.format = 'bioformats';
bf.FrameTimes = 1:Ttotal;
bf.size = [imsize Zdepth];
video = bf;
end

function axeschosen(gcbo,eventdata,selectedChannel)
set(gcf,'UserData',selectedChannel);
uiresume;
end


