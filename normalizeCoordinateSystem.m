function [K,X] = normalizeCoordinateSystem(varargin)
% normalize the coordinate system given by X.
% ------------------------------------------------------------------------
% first check the size of X and fix it if needed
% ------------------------------------------------------------------------
X = varargin{1};
if isnumeric(X)
    x = X;
    y = varargin{2};
    X = [x(:) y(:)];
    ndims=2;
end
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
        X = [x y z]; % build the X matrix
    else % otherwise
        X = [x y]; % make X out of x and y
    end
end

szX = size(X); % check size of X
if any(szX(2) == [2 3])% X is in wront direction ,so transpose it
    X = X'; % transpose it
    szX = size(X); % check size again
end
if szX(1) == ndims % x is not in homogeneous coordinates
    X(ndims+1,:) = ones(1,szX(2)); % add a final row for homogenous coordinates
end


% ------------------------------------------------------------------------
% then find a linear operator K to normalize the cordinate system of X
% ------------------------------------------------------------------------
d = mean(X,2); % find the mean in each direction
if false
    m = [1 1 1];
else
    m = max(X,[],2)-min(X,[],2); % find the dimensions
end
if ndims ==2
    K.inv = [m(1) 0    d(1)   % build the normalizing matrix
        0    m(2) d(2)
        0    0    1   ];
else
    K.inv = [m(1)  0     0     d(1)   % build the normalizing matrix in 3d
        0     m(2)  0     d(2)
        0     0     m(3)  d(3)
        0     0     0     1   ];
end
K.inv = K.inv;
K.fwd = inv(K.inv); % invert it (we'll need this a lot so just do it once and save it!)

% ------------------------------------------------------------------------
% also compute the centered coordinate system (for jacobians)
% ------------------------------------------------------------------------

%X = K.fwd*X;
X = K.fwd*X;
end
