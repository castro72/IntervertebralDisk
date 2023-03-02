function [regions,centers,polygon] = meshRegions(varargin)
% this function takes in an image and options, then allows a user to select
% a region for analysis, along with a region for exclusion (if desired). it
% outputs the initial points in a numberOfPoints x 2 matrix, with  y x
% (image) coordinates.

% process general inputs
im = varargin{1}; % set the image
imsize = size(im); % get it's size
options = varargin{2}; % get the options

%% switch the Meshing Mode
switch options.Meshing.mode
    case 'boxes' % simple box mode
        if nargin == 2 % check if a polygon was given
            figure; imshow(im); % if it was show the image
            hold on % set up for more plots
            % fill the points in according to the polygon and the spacing
            [regions,centers,polygon] = meshBoxes(options.Meshing.spacing,imsize,...
                options.Meshing.boxSize);
        elseif nargin == 3 % if polygon was given
            polygon = varargin{3}; % just use it
            % fill the points in according to the polygon and the spacing
            % using given polygon
            [regions,centers,polygon] = meshBoxes(options.Meshing.spacing,imsize,...
                options.Meshing.boxSize,polygon);
        end
    case 'points'
        if nargin == 3
            centers = varargin{3};
            regions = makeBoxRegions(centers,options.Meshing.boxSize);
        else
            figure; imshow(im); % if it was show the image
            hold on % set up for more plots
            [regions,centers] = meshPoints(options.Meshing.boxSize); % get the regions by selection
        end
    case 'multiline'
        if nargin == 2 % check if a polygon was given
            figure; imshow(im); % if it was show the image
            hold on % set up for more plots
            % fill the points in according to the polygon and the spacing
            [regions,centers,polygon] = meshBoxesLine(options.Meshing.spacing,...
                options.Meshing.boxSize);
        elseif nargin == 3 % if polygon was given
            points = varargin{3}; % just use it
            % fill the points in according to the polygon and the spacing
            % using given polygon
            [regions,centers,polygon] = meshBoxesLine(options.Meshing.spacing,...
                options.Meshing.boxSize,points);
        end
    case 'centerBoxes' % simple box mode
        [regions,centers] = meshCenterBoxes(options.Meshing.spacing,imsize,...
            options.Meshing.boxSize, options.Meshing.AnalysisRegionPrecentage);
    case '25Dsurface'
                regions = find3Dsurfaces(im,options);
end



if ~isempty(options.Meshing.Mask)
   % we have to mask some points out
   isValid = true(1,size(centers,1));
   for i = 1:size(centers,1)
       ind = sub2ind(size(options.Meshing.Mask),centers(i,1),centers(i,2),centers(i,3));
       isValid(i) = options.Meshing.Mask(ind);
   end
   centers = centers(isValid,:);
   regions = regions(isValid);
end
% get the camera calibration
if isa(options.CameraCalibration,'calibrationSession')
    K = calibrationParameterstoMatrix(options.CameraCalibration);
else
    K  = options.CameraCalibration;
end
% normalize the coordinates in the regions
if options.is3d
    for i = 1:numel(regions) % loop over regions
        % get X,Y,Z values
        X = regions(i).x; Y = regions(i).y; Z = regions(i).z;
        [x,y,z] = transform3D(X,Y,Z,@(X) K.fwd*X); % transform them with the calibraiton
        regions(i).x = x; regions(i).y = y; regions(i).z = z; % save the results
    end
else
    for i = 1:numel(regions)
        X = regions(i).x; Y = regions(i).y;
        fcn = @(X) K.fwd*X;
        [x,y] = transform2D(X,Y,fcn); % transform to normalized image coordinates
        %[regions(i).N] = normalizeCoordinateSystem(x,y); % normalize hte new coordinate system
        % perform additional processing
        switch options.Meshing.mode
            case 'multiline'
                N = regions(i).N;
                R = regions(i).R;
                fcn = @(X) N.inv * R * N.fwd * X;
                [x,y] = transform2D(x,y,fcn);
        end
        regions(i).x = x; regions(i).y = y;
    end
end
%close all
end

function [regions,centers,polygon] = meshBoxes(varargin)
spacing = varargin{1}; % get the spacing
imsize = varargin{2}; % get the image size
boxSize = varargin{3}; % get the box size
if nargin == 4 % check if we the user gave a polygon
    polygon = varargin{4};
    centers = getBoxCenters(polygon,spacing,imsize);
%     regions = makeBoxRegions(centers,boxSize); - HLC edit - i don't think
%     this line is necessary because the same line is on 143 
    finished = true;
else
    finished = false;
end
while ~finished
    % ask the user to draw a polygon around the region of interest
    polygon = polymesh('Create a polygon around the region of interest');
    centers = getBoxCenters(polygon,spacing,imsize); % get the centers in the polygon
    bxs = showBoxes(centers,boxSize);
    % removed exclusion for later
    %     h1 = plot(x(in),y(in),'.b');
    %     polyInner = polymesh('Create a polygon around the region to exclude');
    %     yv = polyInner(:,2);
    %     xv = polyInner(:,1);
    %     out = inpolygon(x,y,xv,yv);
    %     h2 = plot(x(out),y(out),'.r');
    %     h3 = plot(x(in & ~out),y(in & ~out),'.g');
    %     points = [y(in & ~out) x(in & ~out) (1:sum(in & ~out))'];
    % nPoints = sum(in & ~out);
    
    choice = questdlg('Accept or repeat?','Confirm Selection','Accept','Repeat','Accept');
    % Handle response
    switch choice
        case 'Accept'
            finished = true;
            close(gcf)
        case 'Repeat'
            finished = false;
            delete(h1); delete(h2); delete(h3);
    end
end
regions = makeBoxRegions(centers,boxSize);
end

function [regions,centers] = meshPoints(boxSize)
finished = false;
while ~finished
    [x,y] = getpts;
    centers = round([y x]);
    regions = makeBoxRegions(centers,boxSize);
    bxs = showBoxes(centers,boxSize);
    choice = questdlg('Accept or repeat?','Confirm Selection','Accept','Repeat','Accept');
    % Handle response
    switch choice
        case 'Accept'
            finished = true;
            close(gcf)
        case 'Repeat'
            finished = false;
            delete(bxs);
    end
end
end

function centers = getBoxCenters(polygon,spacing,imsize)
xv = polygon(:,1); % get the x points of the polygon
yv = polygon(:,2); % get the y points of the polygon
% create a grid of x and y points over the whole image according to
% spacing
x = 1:spacing(2):imsize(2); % create x positions
y = 1:spacing(1):imsize(1); % create y positions
[x,y] = ndgrid(x,y);
x = x(:); y = y(:); % linearize x and y
in = inpolygon(x,y,xv,yv); % check which are within the polygon
centers = [y(in) x(in)]; % arrange the result
end

function bxs = showBoxes(centers,boxSize)
% show boxes
% create a vertices array around 0
dist = [(boxSize(1)-1)/2 (boxSize(2)-1)/2];
v = [-dist(1)  dist(1)  dist(1) -dist(1)
    -dist(2) -dist(2)  dist(2)  dist(2)];
% preallocate arrays
nPoints = size(centers,1); % get the number of points
nv = zeros(4,2,nPoints); % preallocate the vertices
bxs = zeros(1,nPoints);
% draw boxes around each point

for i = 1:nPoints
    nv(:,:,i) = v'+repmat(fliplr(centers(i,:)),4,1);
    bxs(i) = fill(nv(:,1,i),nv(:,2,i),'r','FaceColor','none');
end
end

function poly = polymesh(user_message)
title(user_message)
new_poly = impoly;
poly = wait(new_poly);
%delete(new_poly)
end


function [regions,centers] = meshCenterBoxes(spacing,sz,boxSize,percent)
ndims = numel(spacing);
if ndims == 3
    is3d = true;
else
    is3d = false;
end
center = sz(1:ndims)/2; % find the center of the image
width = sz(1:ndims).*real(percent); % find the width of the region we're going to track

edges = round([center-width/2 ; center+width/2]); % find the box we're going to track
xv = edges(1,1):spacing(1):edges(2,1); % get the x center values
yv = edges(1,2):spacing(2):edges(2,2); % get the y center values
if is3d
    zv = edges(1,3):spacing(3):edges(2,3); % get the z center values
end
% adjust if the edges are less than the range
if range(edges(:,1))<spacing(1)
    xv = round(center(1));
end
if range(edges(:,2))<spacing(2)
    yv = round(center(2));
end
if is3d

    if range(edges(:,3))<spacing(3) && is3d
        zv = round(center(3));
    end

    [x,y,z] = ndgrid(xv,yv,zv); % grid the center vectors
else
    [x,y] = ndgrid(xv,yv);

end
% remove centers from center
if ~isreal(percent) % we have some points to remove from center (has an imaginary part)
    isvalid = false(1,numel(x));
    r = mean((imag(percent).*sz(1:ndims)));

    for i = 1:numel(x)
        
        if is3d
            d = norm([[x(i) y(i) z(i)] - center]);
        else
            d = norm([[x(i) y(i)] - center]);
        end
        isvalid(i) = d>r;
    end
    x(~isvalid) = [];
    y(~isvalid) = [];
    if is3d
        z(~isvalid) = [];
    end
end
if is3d
    centers = [x(:) y(:) z(:) x(:) y(:) z(:)]; % make an array of the centers
    nd = 3;
else
    centers = [x(:) y(:) x(:) y(:)]; % make an array of the centers
    nd = 2;
end

nregions = size(centers,1); % find out how many regions we have
v = [-(boxSize(1:nd)-1)/2 (boxSize(1:nd)-1)/2]; % create a vertices array
vertices = centers+repmat(v,nregions,1); % find the vertices of each region
regions = struct(); % preallocate the regions
for n = 1:nregions % loop over each region
    c = vertices(n,:); % get the specific vertices for that region

    xv = c(1):c(ndims+1); % create x vector of points
    yv = c(2):c(ndims+2); % create y vector of points
    if is3d
        zv = c(3):c(ndims+3); % create z vector of points
        [x,y,z] = ndgrid(xv,yv,zv); % grid the points
        regions(n).x = x; regions(n).y = y; regions(n).z = z; % save the region
    else
        [x,y] = ndgrid(xv,yv); % grid the points
        regions(n).x = x; regions(n).y = y; % save the region
    end
end
centers = centers(:,1:ndims); % truncate the centers for saving

end


function [regions,centers,poly] = meshBoxesLine(spacing,boxSize,varargin)
if nargin < 3
    new_poly = impoly(gca,'Closed',false);
    poly = wait(new_poly);
    close all
else poly = varargin{1};
end
man = manifold(poly,'spacing',spacing(1),'method','pchip');
centers = [man.y ;man.x]';
regions = makeBoxRegions(centers,boxSize);
R = man.rotationMatrices;
t = man.t;
for i = 1:numel(regions)
    regions(i).R = R(:,:,i)';
    regions(i).t = t(i);
end
end

function regions = find3Dsurfaces(im,options)
%% show pair of images and select a region from one of them
plot(im)
poly = polymesh('Select polygon in either image');
sz = size(im);
imSelected = 0;
while imSelected == 0
    if all(poly(:,1)<sz(2))
        imSelected = 1;
        imOther = 2;
    elseif all(poly(:,1)>sz(2))
        imOther = 1;
        imSelected = 2;
         poly(:,1) = poly(:,1)-sz(2);
    else
        warning('Selected pixels from 2 images, please reselect only from left or right');
        poly = polymesh('Select polygon in either image');
    end
end
xv = poly(:,1); % get the x points of the polygon
yv = poly(:,2); % get the y points of the polygon
ROI = [min(xv) min(yv) max(xv)-min(xv) max(yv)-min(yv)];

%% detect features and find a projective transform between the two cameras
pointsSelected = detectSURFFeatures(im{imSelected},'ROI',ROI);
x = pointsSelected.Location(:,1);
y = pointsSelected.Location(:,2);
in = inpolygon(x,y,xv,yv); % check which are within the polygon
pointsSelected = pointsSelected(in);
if length(pointsSelected)>10 % if we have enough points
    pointsOther = detectSURFFeatures(im{imOther});
    [featuresSelected,validSelected] = extractFeatures(im{imSelected},pointsSelected);
    [featuresOther,validOther] = extractFeatures(im{imOther},pointsOther);
    indexPairs = matchFeatures(featuresSelected,featuresOther);
    matchedPointsSelected = validSelected(indexPairs(:, 1), :);
    matchedPointsOther = validOther(indexPairs(:, 2), :);
else % we don't have enough points, lets just toss some points onto the surface and use NXC to find them
    fprintf('Finding point correspondences using nxc... this make take a few minuntes...\n')
    tic
    pointCorrespondences = nxcMeshAndMatch(im,poly,imSelected);
    t = toc;
    fprintf('done! %d elapsed\n',t);
end
% use RANSAC to eliminate outliers between matched points and the surface
[~,ip1,ip2] = estimateGeometricTransform(pointCorrespondences(:,1:2),pointCorrespondences(:,3:4),'projective');

%% find the surface and create a world coordinate system based on it
% find world positions relative to camera 1
worldPoints = triangulate(ip1,ip2,options.CameraCalibration);
% find normal vectors of the plane by fitting the points to a plane
worldCenter = mean(worldPoints);
wp = worldPoints-repmat(worldCenter,size(worldPoints,1),1);
% normal*x = 0
[U,S,V] = svds(wp);
% reproject world points into plane
S(3,3) = 0;
wp = U*S*V';
% project relative world points onto the planar coordinate system of the surface
Xsurface = wp*V;
% compute homographies relating the new world coordinate system to each camera
[homographies,~] = computeHomographies(cat(3,ip1,ip2), Xsurface(:,1:2));
% find camera 1 pose
K(:,:,1)= options.CameraCalibration.CameraParameters1.IntrinsicMatrix';
[rotationVectors(1,:), translationVectors(1,:)] = computeExtrinsics(K(:,:,1), homographies(:,:,1));
% find camera 2 pose
K(:,:,2) = options.CameraCalibration.CameraParameters2.IntrinsicMatrix';
[rotationVectors(2,:), translationVectors(2,:)] = computeExtrinsics(K(:,:,2), homographies(:,:,2));
% estimate pose between the cameras based on homographies
H.ItoII = homographies(:,:,2)/homographies(:,:,1);
H.IItoI = homographies(:,:,1)/homographies(:,:,2);
[rVcam.ItoII,tVcam.ItoII] = computeExtrinsics(K(:,:,1)\K(:,:,2), H.ItoII);
[rVcam.IItoI,tVcam.IItoI] = computeExtrinsics(K(:,:,2)\K(:,:,1), H.IItoI);
%% mesh points on the world surface -- needs some cleanup work
% caluclate the grid based on spacing input
minPoly = min(poly); % find the minimum of the polygon
maxPoly = max(poly); % find the max of the polygon
camCenter = (maxPoly-minPoly)/2+minPoly; % find the cetner of the surface using the average
camGrid = options.Meshing.spacing; % get the spacing input by the user in camera pixels
grid = [[-camGrid/2 ; camGrid/2]'+repmat(camCenter',1,2); 1 1]; % set up a grid of pixels
gridWorld1 = homographies(:,:,1)\grid; gridWorld1= gridWorld1(1:2,:)./repmat(gridWorld1(3,:),2,1); % find the grid locations on the surface in world coordinates
gridWorld2 = homographies(:,:,2)\grid; gridWorld2= gridWorld2(1:2,:)./repmat(gridWorld2(3,:),2,1); % find the grid locations on the surface in world coordinates
worldSpacing = abs(mean([gridWorld1(:,2)-gridWorld1(:,1) gridWorld2(:,2)-gridWorld2(:,1)],2)); % using both cameras determine an appropriate spacing value in the world coordinates
pxToWorld = worldSpacing./camGrid'; % find the conversion factors from pixels to world coordinates
worldBox = options.Meshing.spacing.*pxToWorld'; % create a grid box in the world based on that conversion
worldBox = [-worldBox(1)/2 -worldBox(2)/2 worldBox(1)/2 worldBox(2)/2]; % set up the box as -1/2 + 1/2 the box (box centered about 0)
% push the polygon into the world coordinate system for gridding
polyworld = homographies(:,:,imSelected)\[poly ones(size(poly,1),1)]'; % use selected image 
polyworld = polyworld(1:2,:)./repmat(polyworld(3,:),2,1); % convert from homogenous coordinates
mn = min(polyworld'); % find the min
mx = max(polyworld'); % find the max
xv = mn(1):worldSpacing(1):mx(1); % grid min to max in x
yv = mn(2):worldSpacing(2):mx(2); % grid min to max in y
[x,y] = ndgrid(xv,yv); % create gridded coordinates from grid vectors
in = inpolygon(x,y,polyworld(1,:),polyworld(2,:)); % check which points are within the polygon
x = x(in); % elimate points outside of the polygon
y = y(in); % "" 
X = [x y ones(numel(x),1)]; % set up homogenous coordinate s
grdpts1 = homographies(:,:,1)*X'; % get the gridded points in image 1
centers{1} = grdpts1(1:2,:)./repmat(grdpts1(3,:),2,1); % convert from homogenous coordinates and save as centers
% repeat above for image 2
grdpts2 = homographies(:,:,2)*X';
centers{2} = grdpts2(1:2,:)./repmat(grdpts2(3,:),2,1);
% build each region in the world as [minx miny maxx maxy]
worldBoxes = [x y x y]+repmat(worldBox,size(X,1),1);
% create an lk25region for each resulting region
regions= lk25Region(worldBoxes,pxToWorld,homographies);



end

function pointCorrespondences = nxcMeshAndMatch(im,polygon,selected,options)
% grid the first image
other = setdiff([1 2],selected);
sz = size(im);
sp = [100 100];
xv = [fliplr(round(sz(1)/2):-sp(1):0) (round(sz(1)/2)+10):sp(1):sz(1)];
yv = [fliplr(round(sz(2)/2):-sp(2):0) (round(sz(2)/2)+10):sp(2):sz(2)];
[x,y] = ndgrid(xv,yv); x = x(:); y = y(:);

% remove points outside of polygon
xv = polygon(:,1);
yv = polygon(:,2);
in = inpolygon(x,y,xv,yv); % check which are within the polygon
x = x(in); y = y(in);
% setup bounds
boxSize = [20 20];
bx = [-boxSize(1)/2 -boxSize(2)/2 boxSize(1)/2 boxSize(2)/2];
X = [x y x y];
bxs = X+repmat(bx,size(X,1),1);
% find corresponeces using nxc
temp = im{other};
for i = 1:size(bxs,1)
    test = im{selected}(bxs(i,2):bxs(i,4),bxs(i,1):bxs(i,3));
    C = normxcorr2(test,temp);
    peak = findpeak(C,0.01);
    peaks(i,:) = [peak.y,peak.x]-boxSize/2;
end
if selected == 1
    pointCorrespondences = [x y peaks];
else
    pointCorrespondences = [peaks x y];
end
end

function h = draw25DBoxes(boxVertices)
hold on
clear vert
v = [1 2; 1 4 ; 3 4 ; 3 2];
f = [1 2 3 4];
for i = 1:size(boxVertices,1)
    cVert = boxVertices(i,:);
    vert(:,:,i) = cVert(v);
end
verto = vert;
sz = size(vert);
vert = reshape(permute(vert,[2 1 3]),2,[])';
vert = [vert ones(size(vert,1),1)]';

vertCam1 = homographies(:,:,1)*vert;
vertCam1 = (vertCam1(1:2,:)./repmat(vertCam1(3,:),2,1));
vertCam1 = permute(reshape(vertCam1,sz(2),sz(1),sz(3)),[2 1 3]);

vertCam2 = homographies(:,:,2)*vert;
vertCam2 = (vertCam2(1:2,:)./repmat(vertCam2(3,:),2,1));
vertCam2 = permute(reshape(vertCam2,sz(2),sz(1),sz(3)),[2 1 3]);

vert = vert(1:2,:)';
vert =reshape(vert,sz(1),sz(2),sz(3));
imsize = size(im);
for i = 1:size(vertCam1,3)
    patch('Faces',f,'Vertices',vertCam1(:,:,i),'FaceColor','none','EdgeColor','b');
    v2 = vertCam2(:,:,i);
    v2(:,1) = v2(:,1)+imsize(2);
    patch('Faces',f,'Vertices',v2,'FaceColor','none','EdgeColor','b')
end
    
end