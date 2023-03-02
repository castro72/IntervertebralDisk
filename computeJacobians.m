function dWdP = computeJacobians(x,fcn,p)
debug = false; % turn on to plot the jacobians
userwarp = false; % not needed at the moment, but retained for future warp types
% check if we need to pass a p value (or if we've already computed dWdP via user warp)
if nargin == 2 && ~userwarp
    dWdP = feval(fcn,x); % dWdP is returned in the form (data,p,x)
else
    dWdP = feval(fcn,x,p); % dWdP is returned in the form (data,p,x)
end
if debug % if debugging is enabled show the jacobians in a plot
    showJacobians(dWdP,x); % show them
end
end

function showJacobians(dWdP,X)
% create a regular grid to display the points on
% if X and Y dimensions are small (WILL be the case in small coordinates)
% then we need to approximate how large and small the two dimensions are to 
% get how much to give to each one.
nP = size(dWdP); % get the size
mins = min(X,[],2); % find the mins
maxs = max(X,[],2); % find the maxes
dims = maxs-mins; dims = dims(1:2); % compute the dimensions of the edges
normd = dims/sum(dims); % normalize it
nelements = size(X,3); % find the number of elements we want (based on the number of points)
sides = ceil(sqrt(nelements/prod(normd))*normd); % the number of elements 
% on each edge to use to make a rectangular interpolation region
[x,y] = meshgrid(mins(1):maxs(1),mins(2):maxs(2)); % mesh the grid
ind = 0; % start the counter
for dim = 1:nP(3) % loop over the dimensions
    for P = 1:nP(2) % loop over the parameters
        ind = ind+1; % increment the counter
        subplot(nP(3),nP(2),ind) % get to the right subplot
        % create a scattered interpolant
        F = scatteredInterpolant(X(1,:)',X(2,:)',dWdP(:,P,dim));
        imagesc(F(x,y)) % interpolate and display the jacobian
    end
end
end