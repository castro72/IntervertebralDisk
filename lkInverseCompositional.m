function [mechanics,preComp,regions,nxc] = lkInverseCompositional(video,options) 
%% generalized LK code - no requirements on the type of warp or the
% structure of the regions
% the warp will be specified as a function in options.lk.warp
%
% Inputs Description:
% video - a descriptor indentifying the video we are using. will be
%   interpreted by getVideoFrame function
% options - a structure containing all options related to image processing,
%   optimization, and outputs.

%% process the video
% -------------------------------------
% load the first frame with gradients
% -------------------------------------
getGradients = true;
[template,RGBImage,gradients,options] = getVideoFrame(video,video.Frames(1),options,getGradients);
size(RGBImage)
% -------------------------------------
% mesh the template image
% -------------------------------------
[regions,X] = meshRegions(RGBImage,options);

% -------------------------------------
% perform optional normal cross correlation (for comparison only!)
% -------------------------------------
tic
if options.nxc 
    nxc = nxcSolver(template,regions,video,options);
else nxc = [];
end
nxc.nxcTime = toc;

% -------------------------------------
% perform precomputation for each region 
% -------------------------------------
tic
preComp = preComputation(gradients,regions,template,options);
pcTime = toc;
tic
% -------------------------------------
% register every region in each frame 
% -------------------------------------
[P,stats] = LKSolver(video,preComp,options);

% -------------------------------------
% process the results to get the mechanics
% -------------------------------------
mechanics = processResults(preComp,video,P,X,options);
mechanics.pcTime = pcTime;
mechanics.LKlooptime = toc;
mechanics.options = options;
mechanics.stats = stats;
end

