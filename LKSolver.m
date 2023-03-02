function [P,stats] = LKSolver(video,preComp,options)
% main LK solution loop over all frames
P = initializeResults(preComp,video.Frames);
% go through each region and perform precomputation, saving the result
if numel(preComp)*prod(options.Meshing.boxSize)/1000^2>10 && options.statusBars; % will take longer than ~5 seconds
    options.waitBar.type = 'regions';
    options.waitBar.h = waitbar(0,['1 /' num2str(numel(preComp))],'Name',['Processing ' num2str(numel(video.Frames)) ...
        ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
elseif options.statusBars
    options.waitBar.type = 'frames';
    options.waitBar.h = waitbar(0,['1 /' num2str(numel(video.Frames))],'Name',['Processing ' num2str(numel(video.Frames)) ...
    ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
options.NumberOfFrames = video.NumberOfFrames;
% start a loop
for f = 2:1:numel(video.Frames) % dont use first frame
    if options.statusBars && ~isempty(getappdata(options.waitBar.h,'canceling')) % check if cancel was pushed on the waitbar
        break % if it was cancel the code execution
    end
    warping = getVideoFrame(video,video.Frames(f),options); % load the next image
    PPrev = P(:,:,f-1); % grab the previous results
    % register each point
    [P(:,:,f),stats(f)] = registerRegions(warping,preComp,PPrev,f,options);
    if any(~stats(f).converged) && options.Optimization.recovery
         sum(~stats(f).converged)
         recoverRegions(stats,video,warping,preComp,PPrev,f,options)
    end
    % if plotlive
    % set(h1,'CData',warping.Values);
    %   plotframe(PC,r,mechanics)
    % end
    % update the waitbar
    if options.statusBars
        switch options.waitBar.type
            case 'frames'
                waitbar(f/numel(video.Frames),options.waitBar.h, ...
                    ['Frames: ' num2str(f) '/' num2str(numel(video.Frames))])
            case 'regions'
                waitbar(0,options.waitBar.h,...
                    ['Frames: ' num2str(f) '/' num2str(numel(video.Frames)) ...
                    '    Regions: 0/' num2str(numel(preComp))])
        end
    end
end
if options.statusBars
    delete(options.waitBar.h)
end
end

function [P,stats] = initializeResults(preComp,frames) % initialize parameter results matrix P
    nRegions = numel(preComp);
    nP = size(preComp(1).VTdWdP,2);
    nFrames = numel(frames);
    P = zeros(nRegions,nP,nFrames);
end