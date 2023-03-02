function [r,stats] = registerRegions(warpImage,PC,rPrev,frame,options)

debug = false; % turn this on to plot each registration
r = rPrev; % preallocate r 
warpingFcn = options.fwdWarpingFcn; % grab the warping function for convenience
composeP = options.composition; % grab the composition function for convenience
hTemp = []; hWarp = []; hDiff = []; %#ok<NASGU> % set a null hTemp for debugging
collectStats = nargout == 2; % if the program is asking for 2 outputs, also collect stats on the optimization
if collectStats % stats are requested, preallocate the arrays for them
    dPAll = r; % should be same size as r
    iterAll = zeros(1,numel(PC)); % just an array of numbers
    imerrorAll = zeros(1,numel(PC)); % just an array of numbers
    converged = ones(1,numel(PC)); % preallocate an array for convergence
    status(numel(PC),1) = checkOptimzation(1,zeros(size(rPrev(1,:))),zeros(size(rPrev(1,:))),options);
end
is3d = options.is3d;
%% loop through each point and process it
loopTime = tic; % calculate the time it takes to loop
if options.statusBars
    wBarRegions = strcmp(options.waitBar.type,'regions');
end
for i = 1:numel(PC) % loop through each point    
    Hinv = PC(i).Hinv; % get the inverse Hessian
    if options.Optimization.LM
        H = PC(i).H;  % get the inverse Hessian
        HLM = PC(i).HLM; % get the LM Hessian
        delta = 0.01; % set delta = 0.01 for LM calculation.
    end
    VTdWdP = PC(i).VTdWdP; % get VTdWdP
    P = rPrev(i,:); % grab the previous parameters
    X = PC(i).X.x; % get the original X coordinates
    Y = PC(i).X.y; % get the original Y coordinates
    if is3d % if we have 3d data
        Z = PC(i).X.z; % also get the z coordinate
    end
    temp = PC(i).template; % get template
    iterCount = 0; % initialize an iteration counter
    stop = false; % set stop to off so that the opimization initializes
    failed = false; % set failed to false so that the opimization initializes
%      keyboard;
    if debug % if debugging initialize the debug plot
      [hTemp,hWarp,hDiff] = initializeDebugPlot(temp,hTemp,hWarp,hDiff); %#ok<UNRCH>
    end
    while ~stop && ~failed
        % increment the iteration counter
        iterCount = iterCount+1;
        P_ = P; % save the good update of P
        if is3d
            % Transform X and Y using the warping function with parameters P
            % and normalization function N
            [x,y,z] = transform3D(X, Y, Z, warpingFcn, P , PC(i).N);
            % then take those x y coordinates and grab what they are from the
            % warped image.
            warped = warpImage(x, y, z);
        else
            [x,y] = transform2D(X, Y, warpingFcn, P , PC(i).N);          
            % then take those x y coordinates and grab what they are from the
            % warped image.
            warped = warpImage(x, y);
        end        
        imerror = warped-temp; % calculate the image error
        if options.Optimization.LM % Levenberg Marquardt Optimization (currently not working)
            epsstar = sum(abs(imerror(:)));
            if iterCount == 1 % first run
                eps = epsstar; % save the image error
                % calculate a parameter update
                dP =  (H + delta*HLM) \ VTdWdP' * imerror(:);
                % compose the warps as normal
                P = feval(composeP,P_,options.Optimization.phi * dP);
                % check optimization
                [stop,failed] = checkOptimzation(iterCount,dP,options);
            elseif eps < epsstar % image erorr increased
%                 keyboard
                % fprintf('error increased\n');
                delta = delta*10; % increase delta
                % calculate a new dP
                dP =  (H + delta*HLM) \ VTdWdP' * imerror(:);
                % update the parameters P using the last good update P_
                P = feval(composeP,P_,dP);
                % check if we're stopped (but force dP to be large)
                [stop,failed] = checkOptimzation(iterCount,ones(size(dP)),options);
            elseif eps > epsstar
                delta = delta/10; % decrease delta
                dP =  (H + delta*HLM) \ VTdWdP' * imerror(:); % calculate a new update
                % update the parameters P using the last good update P_
                P = feval(composeP,P_,dP);
                P_ = P; %#ok<NASGU> % save the good update
                eps = epsstar; % save the current image error.
                [stop,failed] = checkOptimzation(iterCount,dP,options);
            end            
        else % standard Guass Newton optimization
            % calculate a parameter update:
            %      dP = H^-1 * VTdWdP^T * (warped - template)
        	dP =  Hinv * VTdWdP' * imerror(:);
            % produce a set of new parameters by composing P and dP -> P
            P = feval(composeP,P,dP);
            % check if the optimization is complete
            status(i) = checkOptimzation(iterCount,P,dP,options);
            stop = status(i).converged;
            failed = status(i).failed;
        end
%         keyboard;
        if debug % if debugging is enabled, update the plots
            updateDebugPlots(warped,imerror,hWarp,hDiff); %#ok<UNRCH>
            pause(.5);
        end   
    end
    % collect results for point
    if failed && collectStats 
        converged(i) = false;
%          keyboard;
    elseif iterCount == 1 % if only one iteration passed, consider the regions identical
        r(i,:) = P; % and use the previous values of P
    else % otherwise use the updated value
        r(i,:)= P;
    end
    % gather optimization statistics for point
    if nargout ==2
        dPAll(i,:) = dP';
        iterAll(i) = iterCount;
        imerrorAll(i) = sum(abs(imerror(:)));
    end
    if options.statusBars && wBarRegions
        if ~isempty(getappdata(options.waitBar.h,'canceling')) % check if cancel was pushed on the waitbar
           break % if it was cancel the code execution
        end
    end
    if options.statusBars && rem(i,30)==0 && wBarRegions 
        waitbar(i/numel(PC),options.waitBar.h,['Frames: ' num2str(frame) '/' num2str(options.NumberOfFrames) ...
                 '    Regions: ' num2str(i) '/' num2str(numel(PC))])
    end
end

%% cleaning up and displaying the results
frameTime = toc(loopTime); % save the frame time
if nargout == 1 % check how many output arguments we have, if there's only one output a simple result
    fprintf('Frame : %d \t Processed %d regions in %0.3f seconds\n',frame,numel(PC),frameTime);
else % otherwise collect optimization statistics
    epsilonMean = mean(sum(abs(dPAll),2)); % get the mean of the epsilon values
    epsilonMin = min(sum(abs(dPAll),2));
    iterMean = mean(iterAll); % get the mean number of iterations
    iterMax = max(iterAll); % get hte max number of iterations
    meanImError = mean(imerrorAll); % find the mean image error
    %fprintf('Frame: %d \t Elapsed: %0.3f \t Epsilon: %f \t Iterations: %0.1f \t Max: %d \t ImageError: %f\n',frame,frameTime,epsilonMean,iterMean,iterMax,meanImError);
    fprintf('Frame: %d \t Elapsed: %0.3f \t Epsilon: %0.4f \t Min: %0.4f \t Iterations: %0.1f \t Max: %d \t ImageError: %0.0f\n',frame,frameTime,epsilonMean,epsilonMin,iterMean,iterMax,meanImError);    
    stats.dP = dPAll; % save the dP values
    stats.iterations = iterAll; % save the iterations
    stats.imageError = imerrorAll; % save the imageerrors
    stats.converged = converged; % save if the region converged
    stats.LoopTime  = frameTime;
    stats.status = status;
end
end

function status = checkOptimzation(iterCount,P,dP,options)
o = options.Optimization; % optimoptions
% -------------------------------------
% check if the optimization is complete
% -------------------------------------
convergence = norm(dP) < o.epsilon; % check if dP < epsilon

% -------------------------------------
% check if the opimization has failed
% -------------------------------------]
if options.is3d
    nd = 3;
else nd = 2;
end
maxIters = iterCount > o.maximumIterations; % check if we're above maximum interations
W = options.warpingFcn(P);

mCondValue = rcond(W(1:nd,1:nd));
mCondStatus = o.minMatrixCondition > mCondValue;

% -------------------------------------
% compile the results
% -------------------------------------
status.converged = convergence;
status.failed = maxIters; % || mCondStatus;
status.mCond = mCondValue;
status.mCondStatus = mCondStatus;
end

function [hTemp,hWarp,hDiff] = initializeDebugPlot(temp,hTemp,hWarp,hDiff)
if isempty(hTemp) % set up a loop to display the image registration
    subplot(1,3,1); hTemp = imagesc(temp(:,:,1));
%     subplot(1,3,1); hTemp = imagesc(temp); % use temp as a default to
%     initialize    - orig
%     subplot(1,3,2); hWarp = imagesc(temp); % ""    - orig
    subplot(1,3,2); hWarp = imagesc(temp(:,:,1)); % "" 
%     subplot(1,3,3); hDiff = imagesc(temp); % ""    - orig
    subplot(1,3,3); hDiff = imagesc(temp(:,:,1)); % ""
else
%   set(hTemp,'CData',temp)    - orig
    set(hTemp,'CData',temp(:,:,1)) % if the loop is already active just update the template image
end
end

function updateDebugPlots(warped,imerror,hWarp,hDiff)
set(hWarp,'CData',warped)
set(hDiff,'CData',imerror)
drawnow;
pause(0.01);
end

