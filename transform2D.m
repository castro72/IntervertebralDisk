function [x,y] = transform2D(X,Y,fct,varargin)

% check inputs
szX = size(X);
szY = size(Y);
if ~all(szX==szY)
    error('X and Y are different sizes') 
end
sz = szX; % set size just equal to x size

%% process and transform the coordinates
% set up
X = X(:); % linearize X
Y = Y(:); % linearize Y
l = ones(numel(X),1);
X = [X(:) Y(:) l]';

% evaluate
if isempty(varargin)
    x_ = feval(fct,X);
else
    x_ = feval(fct,X,varargin);
end

% reshape the results
x = reshape(x_(1,:),sz(1),sz(2));
y = reshape(x_(2,:),sz(1),sz(2));
