function PC = preComputation(gradients,regions,template,options)
%% precomputation function 
% goes through each region and computes everything 
% check if regions are of consistent size
consistentDimensions = options.Meshing.ConsistentDimensions; % check if consistent dimensions
warpType = options.WarpType; % get the warp type
numberOfRegions = numel(regions); % check number of regions

if consistentDimensions % if we have consistent dimensions compute a single Jacobian
    x = regions(1); % get the first region
    [~,x] = normalizeCoordinateSystem(x); % center it's coordinates x
    % compute it's jacobian. use this for all other precomputations
    dWdP = computeJacobians(x,options.jacobianFcn); 
end
nRstr = num2str(numberOfRegions); % going to need this string a lot, just save it
% go through each region and perform precomputation, saving the result
if options.statusBars
    h = waitbar(0,['1 /' nRstr],'Name',['Precomputing ' nRstr ...
        ' regions...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
for idx = 1:numberOfRegions
    % Check for Cancel button press
    if options.statusBars && ~isempty(getappdata(h,'canceling'))
        break
    end
    if consistentDimensions % check if all regions are consistent sizes
        PC(idx) = preComputeRegion(gradients,template,regions(idx),options.is3d,dWdP); % if they are the jacobian was calculated above, pass it along
    else % if regions are not same size, jacobians have to be calculated for each region, don't pass a jacobian
        PC(idx) = preComputeRegion(gradients,template,regions(idx),options.is3d,options.jacobianFcn);
    end
    if options.statusBars
        waitbar(idx/numberOfRegions,h,[num2str(idx) '/' nRstr]) % update waitbar
    end
end
if options.statusBars
    delete(h)
end
end

function r = preComputeRegion(gradients,template,X,is3d,jacobian)
%%
% ------------------------------------------------------------------------
% check number of inputs to determine if we need to calculate a jacobian for this region
% ------------------------------------------------------------------------
if isa(jacobian,'function_handle') % jacobians are not specified, so we need to calculate them
    computeAllJacobians = true;
else 
    computeAllJacobians = false;
    dWdP = jacobian;
end
% ------------------------------------------------------------------------
% center and compute coordinate system of X
% ------------------------------------------------------------------------
if isfield(X,'N')
    if isstruct(X) % check the X input
        x = X.x; x = x(:); % grab the x and y coordinate s
        y = X.y; y = y(:);
        if isfield(X,'z') % check if we have a z variable
            ndims = 3; % if so set ndims to 3
        else % otherwise set it to 2
            ndims = 2;
        end
        if ndims == 3 % if z is a coordinate also include it
            z = X.z; z = z(:); % add z
            Xi = [x y z]; % build the X matrix
        else % otherwise
            Xi = [x y]; % make X out of x and y
        end
    end
    
    szX = size(Xi); % check size of X
    if any(szX(2) == [2 3])% X is in wront direction ,so transpose it
        Xi = Xi'; % transpose it
        szX = size(Xi); % check size again
    end
    if szX(1) == ndims % x is not in homogeneous coordinates
        Xi(ndims+1,:) = ones(1,szX(2)); % add a final row for homogenous coordinates
    end
    szX = size(Xi); % check size of X
    if any(szX(2) == [2 3])% X is in wront direction ,so transpose it
        Xi = Xi'; % transpose it
        szX = size(X); % check size again
    end
    if szX(1) == ndims % x is not in homogeneous coordinates
        Xi(ndims+1,:) = ones(1,szX(2)); % add a final row for homogenous coordinates
    end
    N = X.N;
    x = N.fwd*Xi;
else
    [N,x] = normalizeCoordinateSystem(X); 
end
% calculate image coordinates

% ------------------------------------------------------------------------
% if we need to compute a jacobian for the region compute it
% ------------------------------------------------------------------------
if computeAllJacobians % if jacobians are not given they must be computed!
    dWdP = computeJacobians(x,jacobian); % NOTE will need to UPDATE this if we recompute a region 
end

% ------------------------------------------------------------------------
% perform precomputation
% ------------------------------------------------------------------------
% get the gradient pixels corresponding to the region (image coordinate
% system)

if is3d
    VTx = gradients.x(X.x,X.y,X.z); VTx = VTx(:); 
    VTy = gradients.y(X.x,X.y,X.z); VTy = VTy(:);
    VTz = gradients.z(X.x,X.y,X.z); VTz = VTz(:);
else
    VTx = gradients.x(X.x,X.y); VTx = VTx(:); 
    VTy = gradients.y(X.x,X.y); VTy = VTy(:);
end
% compute del gradient * jacobian (VTdWdP)
xn = numel(VTx); % get number of x elements
nP = size(dWdP,2); % get number of parameters
VTdWdP = zeros(xn,nP);
if is3d
    for j = 1:xn % VT*dW/dp in equations 35 and 36
        VTdWdP(j,:) = [VTx(j)        VTy(j)        VTz(j)     ] *  ...
                      [dWdP(j,:,1) ; dWdP(j,:,2) ; dWdP(j,:,3)];
    end
else
    for j = 1:xn % VT*dW/dp in equations 35 and 36
        VTdWdP(j,:) = [VTx(j) VTy(j)] * [dWdP(j,:,1) ; dWdP(j,:,2)];
    end
end
% compute hessian
H = VTdWdP'*VTdWdP;
HLM = diag(sum(VTdWdP.^2));
% compute inverse hessian
Hinv = inv(H);
% grab the template image
if is3d
    r.template = template(X.x,X.y,X.z);
else
    r.template = template(X.x,X.y);
end
%% save outputs
r.N = N; % save conversion matrix for the specfic region
% r.dWdP = dWdP; % save jacobian ***DISABLED FOR NOW***
r.VTdWdP = VTdWdP; % save del template * jacobian
r.H = H; % save hessian matrix
r.HLM = HLM; % save Levenberg-Marquadt Hessian matrix
r.Hinv = Hinv; % save inverted Hessian
r.X.x = X.x; % save the x coordinates
r.X.y = X.y; % save the y coordinates
if is3d
    r.X.z = X.z;
end
if isfield(X,'R')
    r.R = X.R; % save rotation matrix too
end
if isfield(X,'t')
    r.t = X.t; % save parameterized distance t
end
end