function options = warpingFunctions(options)
% this function contains all of the known warping functions to be used. if
% a new type is needed please add the warping function here, the forward
% warping function, and the compositional function 
%
% specify the warping functions as follows:
%                     
% options.fwdWarpingFcn: warping function W(x;p) found in main equation
% options.composition  : how to compose W(x;p+dp) with differential warping
%                        parameters dp
% options.jacobianFcn  : funciton describing how to compute dW/dP
% options.warpingFcn   : describe the basic linear function (used to simplify 
%                        fwd and composition in this function)
%                        for affine deformations, this is also later used
%                        to do various mechanics analyses

switch options.WarpType
    case {'affine2d' 'affine1d' 'affine1dmanifold'}
        options.warpingFcn = @(P) affine2d(P);
        options.fwdWarpingFcn = @(P,x) affine2dFwd(P,x);
        options.composition = @(P,dP) affine2dComposition(P,dP);
        options.jacobianFcn = @affine2dJacobian;
    case 'affine3d'
        options.warpingFcn = @(P) affine3d(P);
        options.fwdWarpingFcn = @(P,x) affine3dFwd(P,x);
        options.composition = @(P,dP) affine3dComposition(P,dP);
        options.jacobianFcn = @affine3dJacobian;
    case 'affine3dPrinciple'
        options.warpingFcn = @(P) affine3dPrinciple(P);
        options.fwdWarpingFcn = @(P,x) affine3dFwdPrinciple(P,x);
        options.composition = @(P,dP) affine3dCompositionPrinciple(P,dP);
        options.jacobianFcn = @affine3dJacobianPrinciple;
        
    % -------------------------
    % undeveloped functions
    % -------------------------
    case 'projective2d' % projective 2D transformation 
        options.jacobianFcn = @projective2dJacobian; % pass the appropriate function
        options.warpingFcn = @(P) projective2d(P);
        options.fwdWarpingFcn = @(P,x) projective2dFwd(P,x);
        options.composition = @(P,dP) projective2dComposition(P,dP);
    case 'frt' % warp consisting of a deformation (F), a rotation about x,y, and z (R)
               % and a translation in x, y, and z.
        options.jacobian = @frt; % pass the appropriate warp
    otherwise % if we don't have one of these
        if isa(warp,'struct') % check if we have some sort of warp structure
%            userwarp = true; % if we do we're using a user defined warping function
            switch warp.dims % see how many dimensions the user defined function has
                case 2 % if 2 lets compute a 2D user warp
                    if nargin == 2 % p is not given so 
                        dWdP = userWarp2D(x,warp.syms); %#ok<NASGU>
                    else 
                        dWdP = userWarp2D(x,warp.syms,p); %#ok<NASGU> % p is given so pass it
                    end
            end
        else % this isn't a given option for a warp, throw an angry error
            error('no known jacobian or user given warp!')
        end
end
end

%% Affine 2d functions
function W = affine2d(P)
% affine 2d warping function
W = [1+P(1)     P(3)    P(5)
       P(2)   1+P(4)    P(6)
       0        0       1];
end

function x = affine2dFwd(X,inputs)
narginchk(2,2) % check to make sure 2 inputs are passed
P = inputs{1}; % set P = to the first argin
if numel(inputs)==1
    K = eye(3);
else
    K = inputs{2};
end
W = affine2d(P);
x = K.inv*W*K.fwd*X;
end

function P = affine2dComposition(P,dP)
% This composition is the same as building a warp for P and dP:
% W(x,P):                           W(x,dP):
% [1+P(1)  P(3)    P(5) ;           [1+dP(1)  dP(3)    dP(5) ;
%  P(2)    1+P(4)  P(6) ;            dP(2)    1+dP(4)  dP(6) ;
%  0       0       1   ];            0       0       1      ];
% from the parameters P and dP, then muliplying the resulting matrices:
% W(x,P)*W(x,dP)^-1
%
% derived from equation 16
%
%Here it is calculated explicity for speed:

% den = ( (1+dP(1))*(1+dP(4))- dP(2)*dP(3) ) ^-1;
% 
% dP = [-dP(1) - dP(1)*dP(4) + dP(2)*dP(3)
%                     -dP(2)
%                     -dP(3)
%       -dP(4) - dP(1)*dP(4) + dP(2)*dP(3)
%       -dP(5) - dP(4)*dP(5) + dP(3)*dP(6)
%       -dP(6) - dP(1)*dP(6) + dP(2)*dP(5)];
%   
% dP = den*dP;
%   
% P =[P(1) + dP(1) + P(1)*dP(1) + P(3)*dP(2) ...
%     P(2) + dP(2) + P(2)*dP(1) + P(4)*dP(2) ...
%     P(3) + dP(3) + P(1)*dP(3) + P(3)*dP(4) ...
%     P(4) + dP(4) + P(2)*dP(3) + P(4)*dP(4) ...
%     P(5) + dP(5) + P(1)*dP(5) + P(3)*dP(6) ...
%     P(6) + dP(6) + P(2)*dP(5) + P(4)*dP(6)  ];

PW =  affine2d(P);
dPW = affine2d(dP);
W = PW/dPW;
W = W-eye(3);
W(end,:) = [];
P = W(:)';
end

function dWdP = affine2dJacobian(X,~) % build simple affine 2D jacobians
sz = size(X,2);
O = zeros(sz,1); % zero fill for jacobian
l =  ones(sz,1); % one fill for jacobian
x = X(1,:)'; x = x(:); % extract and linearize x
y = X(2,:)'; y = y(:); % extract and linearize y
% these jacobians are known so we just fill in the data for x and y
dWdP(:,:,1) = [ x O  y O  l O ]; % build jacobian for x
dWdP(:,:,2) = [ O x  O y  O l ]; % build jacobian for y
end

%% Affine 3D functions
function W = affine3d(P)
% affine 3d warping function
W = [1+P(1)     P(4)      P(7)  P(10)
       P(2)   1+P(5)      P(8)  P(11)
       P(3)     P(6)    1+P(9)  P(12)
         0        0       0        1];
end

function x = affine3dFwd(X,inputs)
% foward warping function with camera calibration
% x = K^-1 * W * K * X
narginchk(2,2) % check to make sure 2 inputs are passed
P = inputs{1}; % set P = to the first argin
if numel(inputs)==1
    K.fwd = eye(4);
    K.inv = eye(4);
else
    K = inputs{2};
end
W = affine3d(P);
x = K.inv*W*K.fwd*X;
end

function P = affine3dComposition(P,dP)
% This composition is the same as building a warp for P and dP:
% W(x,P) = [1+P(1)     P(4)      P(7)  P(10)
%             P(2)   1+P(5)      P(8)  P(11)
%             P(3)     P(6)    1+P(9)  P(12)
%             0        0       0       1    ];
% W(x,dP) = [1+dP(1)     dP(4)      dP(7)  dP(10)
%              dP(2)   1+dP(5)      dP(8)  dP(11)
%              dP(3)     dP(6)    1+dP(9)  dP(12)
%              0         0          0      1    ];
% from the parameters P and dP, then muliplying the resulting matrices:
% W(x,P)*W(x,dP)^-1
%
% derived from equation 16

PW =  affine3d(P);
dPW = affine3d(dP);
W = PW/dPW;
W = W-eye(4);
W(end,:) = [];
P = W(:)';
end

function dWdP = affine3dJacobian(X,~) % build simple affine 2D jacobians
sz = size(X,2);
O = zeros(sz,1); % zero fill for jacobian
l =  ones(sz,1); % one fill for jacobian
x = X(1,:)'; x = x(:); % extract and linearize x
y = X(2,:)'; y = y(:); % extract and linearize y
z = X(3,:)'; z = z(:); % extract and linearize z
% these jacobians are known so we just fill in the data for x and y
dWdP(:,:,1) = [ x O O  y O O  z O O  l O O ]; % build jacobian for x
dWdP(:,:,2) = [ O x O  O y O  O z O  O l O ]; % build jacobian for y
dWdP(:,:,3) = [ O O x  O O y  O O z  O O l ];
end

%% Affine 3D simplified functions
function W = affine3dPrinciple(P)
% affine 3d warping function
W = [1+P(1)      0         0   P(4)
         0   1+P(2)        0   P(5)
         0       0     1+P(3)  P(6)
         0        0       0      1 ];
end

function x = affine3dFwdPrinciple(X,inputs)
% foward warping function with camera calibration
% x = K^-1 * W * K * X
narginchk(2,2) % check to make sure 2 inputs are passed
P = inputs{1}; % set P = to the first argin
if numel(inputs)==1
    K.fwd = eye(4);
    K.inv = eye(4);
else
    K = inputs{2};
end
W = affine3dPrinciple(P);
x = K.inv*W*K.fwd*X;
end

function P = affine3dCompositionPrinciple(P,dP)
% This composition is the same as building a warp for P and dP:
% W(x,P) = [1+P(1)     P(4)      P(7)  P(10)
%             P(2)   1+P(5)      P(8)  P(11)
%             P(3)     P(6)    1+P(9)  P(12)
%             0        0       0       1    ];
% W(x,dP) = [1+dP(1)     dP(4)      dP(7)  dP(10)
%              dP(2)   1+dP(5)      dP(8)  dP(11)
%              dP(3)     dP(6)    1+dP(9)  dP(12)
%              0         0          0      1    ];
% from the parameters P and dP, then muliplying the resulting matrices:
% W(x,P)*W(x,dP)^-1
%
% derived from equation 16

PW =  affine3dPrinciple(P);
dPW = affine3dPrinciple(dP);
W = PW/dPW;
W = W-eye(4);
W(end,:) = [];
P = W([1 5 9:12])';
end

function dWdP = affine3dJacobianPrinciple(X,~) % build simple affine 2D jacobians
sz = size(X,2);
O = zeros(sz,1); % zero fill for jacobian
l =  ones(sz,1); % one fill for jacobian
x = X(1,:)'; x = x(:); % extract and linearize x
y = X(2,:)'; y = y(:); % extract and linearize y
z = X(3,:)'; z = z(:); % extract and linearize z
% these jacobians are known so we just fill in the data for x and y
dWdP(:,:,1) = [ x O O   l O O ]; % build jacobian for x
dWdP(:,:,2) = [ O y O   O l O ]; % build jacobian for y
dWdP(:,:,3) = [ O O z   O O l ];
end

function W = projective2d(P)
% doi:10.1023/B:VISI.0000011205.11775.fdfrom 
% W(x,P):                
% [1+P(1)  P(3)    P(5) ;
%  P(2)    1+P(4)  P(6) ;
%  P(7)    P(8)    1   ]; 
W = [1+P(1)     P(3)    P(5)
       P(2)   1+P(4)    P(6)
       P(7)     P(8)      1];
end

function x = projective2dFwd(X,inputs)

narginchk(2,2) % check to make sure 2 inputs are passed
P = inputs{1}; % set P = to the first argin
if numel(inputs)==1
    K = eye(3);
else
    K = inputs{2};
end
W = projective2d(P);
x = W*K.fwd*X;
x = x./repmat(x(3,:),3,1);
x = K.inv*x;
end

function P = projective2dComposition(P,dP)
% affine 2d composition. equation 96 from
% Baker, S., & Matthews, I. (2004). Lucas-Kanade 20 Years On: A Unifying 
% Framework. International Journal of Computer Vision, 56(3), 221�255.
% doi:10.1023/B:VISI.0000011205.11775.fdfrom 
% W(x,P):                           W(x,dP):
% [1+P(1)  P(3)    P(5) ;           [1+dP(1)  dP(3)    dP(5) ;
%  P(2)    1+P(4)  P(6) ;            dP(2)    1+dP(4)  dP(6) ;
%  P(7)    P(8)    1   ];            dP(7)    dp(8)    1      ];
PW =  projective2d(P);
dPW = projective2d(dP);
W = PW/dPW;
W = W-eye(3);
P = W([1 2 4 5 7 8 3 6]);
% dP = inverseProjective(dP);
% n = ( 1 + P(7)*dP(5) + P(8)*dP(6) );
% c(1) = P(1) + dP(1) + P(1)*dP(1) + P(3)*dP(2) + P(5)*dP(7) - P(7)*dP(5) - P(8)*dP(6) ;
% c(2) = P(2) + dP(2) + P(2)*dP(1) + P(4)*dP(2) + P(6)*dP(7);
% c(3) = P(3) + dP(3) + P(1)*dP(3) + P(3)*dP(4) + P(5)*dP(8);
% c(4) = P(4) + dP(4) + P(2)*dP(3) + P(4)*dP(4) + P(6)*dP(8) - P(7)*dP(5) - P(8)*dP(6) ;
% c(5) = P(5) + dP(5) + P(1)*dP(5) + P(3)*dP(6);
% c(6) = P(6) + dP(6) + P(2)*dP(5) + P(4)*dP(6);
% c(7) = P(7) + dP(7) + P(7)*dP(1) + P(8)*dP(2);
% c(8) = P(8) + dP(8) + P(7)*dP(3) + P(8)*dP(4);
% P = c/n;
end

function P = inverseProjective(p)
d = det(projective2d(p));
n = d * ( (1+p(1))*(1+p(4)) - p(2)*p(3) );
P(1) = 1 + p(1) - p(6)*p(8) - n;
P(2) =    -p(2) + p(6)*p(7);
P(3) =    -p(3) + p(5)*p(8);
P(4) = 1 + p(4) - p(5)*p(7) - n;
P(5) =    -p(5) - p(4)*p(5) + p(3)*p(6); 
P(6) =    -p(6) - p(1)*p(6) + p(2)*p(5);
P(7) =    -p(7) - p(4)*p(7) + p(2)*p(8);
P(8) =    -p(8) - p(1)*p(8) + p(3)*p(7);
P = P/n;
end
%         options.fwdWarpingFcn = @(P,x) affine2dFwd(P,x);
%         options.composition = @(P,dP) affine2dComposition(P,dP);
%% various jacobians for future use
function dWdP = projective2dJacobian(X,p)
% affine 2d jacobian. equation 95 from
% Baker, S., & Matthews, I. (2004). Lucas-Kanade 20 Years On: A Unifying 
% Framework. International Journal of Computer Vision, 56(3), 221�255.
% doi:10.1023/B:VISI.0000011205.11775.fdfrom 
if nargin == 1
    p = zeros(1,8);
end
sz = size(X,2); 
O = zeros(sz,1); % zero fill for jacobian
l =  ones(sz,1); % one fill for jacobian
x = X(1,:)';  % extract and linearize x
y = X(2,:)';  % extract and linearize y
% these jacobians are known so we just fill in the data for x and y
n = 1 + p(7)*x + p(8)*y;
x7 = -x.*( (1+p(1))*x +     p(3)*y + p(5) )./n; % build x7 position
x8 = -y.*( (1+p(1))*x +     p(3)*y + p(5) )./n; % build x8 position
y7 = -x.*(     p(2)*x + (1+p(4))*y + p(6) )./n; % build y7 position 
y8 = -y.*(     p(2)*x + (1+p(4))*y + p(6) )./n; % build y8 position
dWdP(:,:,1) = [ x O  y O  l O  x7 x8 ]; % build jacobian for x
dWdP(:,:,2) = [ O x  O y  O l  y7 y8 ]; % build jacobian for y
n = repmat(n,[1 8 2]); % replicate n so we can divide it across
dWdP = dWdP./n; % divide by n (see formula in paper!)
end

function dWdP = frtJacobian(X,p)
if nargin == 1
    p = zeros(1,8); % build p
    % the rotational components MUST be non-zero for conversion to a rotation matrix to work!
    % so we're going to set them very close to 0.
    p(1:3) = 0.00001; 
end
% frt is a known warp but it's too complicated to compute the jacobian. we're going to built 
% it as a symbolic function then pass it to userwarp to find the jacobian
sx = sym('sx',[1 2]); % symbolic notation for x (sx)
P = sym('P',[1 8]); % 3 Parameters for Rotation, 2 For strain, 3 for Translation -- 8 parameters
R  = rotationVectorToMatrix([P(1) P(2) P(3)]); % convert R to a rotation tensor
F = [P(7)+1  0          0    % build our deformation tensor
     0       P(8)+1     0 
     0       0          1];
W = F*R*[sx(1);sx(2);1] + [P(4);P(5);P(6)]; % build the symbolic warp jacobian
W = W(1:2)/W(3); % project the warp into 2D space
dWdP = userWarp2D(X,W,p); % pass it to the userwarp function 
end

function dWdP = userWarp2D(X,W,p)
if nargin == 2
    nP = numel(symvar(W))-2; % number of variables in the warp is given by 
    % numel(symvar). there's also an x y pair for the coordinates so subtract 2
    p = zeros(1,nP); % build an initial value matrix for P
end
sx = sym('sx',[1 2]); % remake symbolic variables
P = sym('P',[1 numel(p)]); % have to remake this as a symbolic variable I guess
dWdPsym = jacobian(W,P); % compute the symbolic jacobian
dWdP0 = subs(dWdPsym,P,p); % substitute P for given P values
x = X(1,:)'; % get the x vector
y = X(2,:)'; % get the y vector
dWdP = zeros(numel(x),numel(p),2); % preallocate!
fprintf('Manually computing Jacobians! This may take a moment\n');
for i = 1:numel(x) % loop over each x y value pair
     % compute the jacobian by substituting the current xy value pair
     dWdP(i,:,:) = double(subs(dWdP0,sx,[x(i) y(i)]))'; 
end
fprintf('Jacobian computed!\n');
end
   