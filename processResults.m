function mechanics = processResults(PC,video,P,X,options)
% switch which type of analysis we're doing
switch options.WarpType
    case 'affine2d' % affine2d analysis
        mechanics = processAffine(PC,video,P,X,options);
    case 'affine3d'
        mechanics = processAffine(PC,video,P,X,options);
    case 'affine3dPrinciple'
        mechanics = processAffine(PC,video,P,X,options);
    case 'affine1dmanifold'
        mechanics = processManifold2d(PC,video,P,X,options);
    case 'projective2d'
        mechanics = processProjective(PC,video,P,X,options);
    case 'affine1d' % affine2d analysis
        mechanics = processAffine1d(PC,video,P,X,options);
end
end



function mechanics = processAffine(PC,video,P,X,options)
nd = size(X,2); % number of dimensions
% get the calibration parameters from the calibration matrix
K = calibrationParameterstoMatrix(options.CameraCalibration); 
% get the number of elements
nelement = size(P,1);
nframe = size(P,3);
% preallocate F, E, and X, FLSF, ELSF, and averages
F = zeros(nd,nd,nelement,nframe); % DDE deformation gradient tensor
E = zeros(nd,nd,nelement,nframe); % DDE strain matrix
x = zeros(nd,nelement,nframe); % x1 x2 positions
EIv = zeros(nd,nd,nelement,nframe); % DDE eigenvectors
EI = zeros(nd,nd,nelement,nframe); % DDE eigenvalues
FLSF = zeros(nd,nd,nelement,nframe); % LSF deformation gradient tensor
ELSF = zeros(nd,nd,nelement,nframe); % LSF strain matrix
ELSFIv = zeros(nd,nd,nelement,nframe); % LSF strain eigenvectors matrix
ELSFI = zeros(nd,nd,nelement,nframe); % LSF strain eigenvalues matrix
FLSF_average = zeros(nd,nd,nframe); % average deformation (LSF over all points)
ELSF_average = zeros(nd,nd,nframe); % average strain      (LSF over all points)
DeltaIv = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvectors
DeltaI = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvalues
Cond = zeros(nelement,nframe);
LSFCond = zeros(nelement,nframe);
% go through each element and each frame and calculate the deformation, the position, and
% the strain
if options.statusBars
    h = waitbar(0,['1 /' num2str(nframe)],'Name',['Processing Mechanics for ' num2str(nframe) ...
        ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
for f = 1:nframe
    if options.statusBars 
        if  getappdata(h,'canceling')% check if the cancel button has been pushed
            break % if so quit
        end
    end
    % -------------------------
    % ----- DDE mechanics -----
    % -------------------------
    for i = 1:nelement
        Ft = feval(options.warpingFcn,P(i,:,f)); % use the Warping Function to generate a deformation gradient tensor for this Region
        xt = Ft*PC(i).N.fwd*K.fwd*[X(i,:) 1]'; % warp the coordinates of X into x for LSF.
        xt = xt./xt(end);
        xt = K.inv*PC(i).N.inv*xt;
        Ft = Ft(1:nd,1:nd); % extract the deformation gradient Tensor (2d)
        F(:,:,i,f) = Ft; % add the deformation to the
        x(:,i,f) = xt(1:nd)./repmat(xt(nd+1),nd,1); % insert temporary position into position array
        E_ = 0.5*(Ft'*Ft-eye(nd)); % calculate DDE strain
        E(:,:,i,f) = E_; % insert strain into array
        [EIv_, EI_] = eig(E_); % get the eigenvectors and eigenvalues of the strain
        EIv(:,:,i,f) = EIv_; % insert the eigenvectors into array
        EI(:,:,i,f) = EI_; % insert the eigenvalues into array
        Cond(i,f) = rcond(Ft); % get condition matrix for DDE
    end
    % -------------------------
    % ----- LSF mechanics -----
    % -------------------------
    xt = x(:,:,f)'; % Transpose x (to match X)
    % find elements that are close to one and other to do LSF mechanics
    LSF_groups = knnsearch(xt,xt,'K',options.mechanics.LSF_Points);
    for i = 1:nelement
        g = LSF_groups(i,:); % grab the closest points
        if ~checkUniqueDim(X(g,:)) && options.mechanics.fixDim % ensure points dont lie in a single line or plane
            g = fixUniqueDim(xt,xt(i,:),X,options.mechanics.LSF_Points); % if they are, fix them
            %numel(g)
        end
        FLSF_ = LSF_F(X(g,:),xt(g,:)); % calculate the deformation
        ELSF_ = 0.5*(FLSF_'*FLSF_-eye(nd)); % calculate strain based on deformation
        FLSF(:,:,i,f) = FLSF_; % insert temp variable into entire array
        ELSF(:,:,i,f) = ELSF_; % insert temp variable into entire array
        [ELSFIv_,ELSFI_] = eig(ELSF_); % get LSF eigenvalues
        ELSFIv(:,:,i,f) = ELSFIv_; % insert the eigenvectors into array
        ELSFI(:,:,i,f) = ELSFI_; % insert the eigenvalues into array
        LSFCond(i,f) = rcond(FLSF_); % get condition matrix for LSF
    end
    % caclulate average strains based on displacements
    FLSF_average_ = LSF_F(X,xt);
    FLSF_average(:,:,f) = FLSF_average_;
    ELSF_average(:,:,f) = 0.5*(FLSF_average_'*FLSF_average_-eye(nd));
    % update waitbar
    if options.statusBars
        waitbar(f/(nframe+1),h,[num2str(f) '/' num2str(nframe)])
    end
end
% -------------------------
% -------- SIMPLE ---------
% -------------------------
if options.statusBars
    waitbar(nframe/(nframe+1),h,'Calculating Delta (SIMPLE)')
end 
Delta = F-FLSF; % get general delta
for f = 1:nframe % loop through each frame
    for i = 1:nelement % loop through each region
        [DeltaIv_,DeltaI_] = eig(Delta(:,:,i,f)); % get Delta eigenvalues
        DeltaIv(:,:,i,f) = DeltaIv_; % insert the eigenvectors into array
        DeltaI(:,:,i,f) = DeltaI_; % insert the eigenvalues into array
    end
end

% -------------------------
% ---- gather results -----
% -------------------------
if options.statusBars
    waitbar(nframe/(nframe+1),h,'Gathering results....')
end
% general
m.x = x; % x1 x2 positions
m.FrameTimes = video.FrameTimes(video.Frames);
m.P = P; % raw data
% DDE
m.DDE.F = F; % DDE deformation gradient tensor
m.DDE.E = E; % DDE strain matrix
m.DDE.EIv = EIv; % DDE eigenvectors
m.DDE.EI = EI; % DDE eigenvalues
m.DDE.cond = Cond; % DDE matrix condition
% LSF
m.LSF.F = FLSF; % LSF deformation gradient tensor
m.LSF.E = ELSF;  % LSF strain matrix
m.LSF.EIv = ELSFIv;  % LSF strain eigenvectors matrix
m.LSF.EI = ELSFI;  % LSF strain eigenvalues matrix
m.LSF.cond = LSFCond; % LSF matrix condition
% SIMPLE
m.SIMPLE.Delta = Delta; % xy Delta
m.SIMPLE.DeltaIv = DeltaIv; % preallocate SIMPLE eigenvectors
m.SIMPLE.DeltaI = DeltaI; % preallocate SIMPLE eigenvalues
% averages
m.FrameAverage.F = FLSF_average;  % average deformation (LSF over all points)
m.FrameAverage.E = ELSF_average;  % average strain      (LSF over all points)
% save the results
mechanics = m;
% kill waitbar
if options.statusBars
    delete(h)
end
end


function mechanics = processProjective(PC,video,P,X,options)
nd = size(X,2); % number of dimensions
% get the calibration parameters from the calibration matrix
K = calibrationParameterstoMatrix(options.CameraCalibration); 
% get the number of elements
nelement = size(P,1);
nframe = size(P,3);
% preallocate F, E, and X, FLSF, ELSF, and averages
F = zeros(nd,nd,nelement,nframe); % DDE deformation gradient tensor
E = zeros(nd,nd,nelement,nframe); % DDE strain matrix
x = zeros(nd,nelement,nframe); % x1 x2 positions
EIv = zeros(nd,nd,nelement,nframe); % DDE eigenvectors
EI = zeros(nd,nd,nelement,nframe); % DDE eigenvalues
FLSF = zeros(nd,nd,nelement,nframe); % LSF deformation gradient tensor
ELSF = zeros(nd,nd,nelement,nframe); % LSF strain matrix
ELSFIv = zeros(nd,nd,nelement,nframe); % LSF strain eigenvectors matrix
ELSFI = zeros(nd,nd,nelement,nframe); % LSF strain eigenvalues matrix
FLSF_average = zeros(nd,nd,nframe); % average deformation (LSF over all points)
ELSF_average = zeros(nd,nd,nframe); % average strain      (LSF over all points)
DeltaIv = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvectors
DeltaI = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvalues
Cond = zeros(nelement,nframe);
LSFCond = zeros(nelement,nframe);
% go through each element and each frame and calculate the deformation, the position, and
% the strain
if options.statusBars
    h = waitbar(0,['1 /' num2str(nframe)],'Name',['Processing Mechanics for ' num2str(nframe) ...
        ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
for f = 1:nframe
    if options.status && Barsgetappdata(h,'canceling')% check if the cancel button has been pushed
        break % if so quit
    end
    
    
    % -------------------------
    % ----- DDE mechanics -----
    % -------------------------
    for i = 1:nelement
        Ft = feval(options.warpingFcn,P(i,:,f)); % use the Warping Function to generate a deformation gradient tensor for this Region
        xt = Ft*PC(i).N.fwd*K.fwd*[X(i,:) 1]'; % warp the coordinates of X into x for LSF.
        xt = xt./xt(end);
        xt = K.inv*PC(i).N.inv*xt;
        T = eye(3); T([3 6]) = Ft([3 6]); Ft = Ft/T;
        Ft = Ft(1:nd,1:nd); % extract the deformation gradient Tensor (2d)
        F(:,:,i,f) = Ft; % add the deformation to the
        x(:,i,f) = xt(1:nd); % insert temporary position into position array
        E_ = 0.5*(Ft'*Ft-eye(nd)); % calculate DDE strain
        E(:,:,i,f) = E_; % insert strain into array
        [EIv_, EI_] = eig(E_); % get the eigenvectors and eigenvalues of the strain
        EIv(:,:,i,f) = EIv_; % insert the eigenvectors into array
        EI(:,:,i,f) = EI_; % insert the eigenvalues into array
        Cond(i,f) = rcond(Ft); % get condition matrix for DDE
    end
    
    
    
    
    % -------------------------
    % ----- LSF mechanics -----
    % -------------------------
    xt = x(:,:,f)'; % Transpose x (to match X)
    % find elements that are close to one and other to do LSF mechanics
    LSF_groups = knnsearch(xt,xt,'K',options.mechanics.LSF_Points);
    for i = 1:nelement
        g = LSF_groups(i,:); % grab the 5 closest points
        FLSF_ = LSF_F(X(g,:),xt(g,:)); % calculate the deformation
        ELSF_ = 0.5*(FLSF_'*FLSF_-eye(nd)); % calculate strain based on deformation
        FLSF(:,:,i,f) = FLSF_; % insert temp variable into entire array
        ELSF(:,:,i,f) = ELSF_; % insert temp variable into entire array
        [ELSFIv_,ELSFI_] = eig(ELSF_); % get LSF eigenvalues
        ELSFIv(:,:,i,f) = ELSFIv_; % insert the eigenvectors into array
        ELSFI(:,:,i,f) = ELSFI_; % insert the eigenvalues into array
        LSFCond(i,f) = rcond(FLSF_); % get condition matrix for LSF
    end
    % caclulate average strains based on displacements
    FLSF_average_ = LSF_F(X,xt);
    FLSF_average(:,:,f) = FLSF_average_;
    ELSF_average(:,:,f) = 0.5*(FLSF_average_'*FLSF_average_-eye(nd));
    % update waitbar
    waitbar(f/(nframe+1),h,[num2str(f) '/' num2str(nframe)])
end
% -------------------------
% -------- SIMPLE ---------
% -------------------------
if options.statusBars 
    waitbar(nframe/(nframe+1),h,'Calculating Delta (SIMPLE)')
end
Delta = F-FLSF; % get general delta
for f = 1:nframe % loop through each frame
    for i = 1:nelement % loop through each region
        [DeltaIv_,DeltaI_] = eig(Delta(:,:,i,f)); % get Delta eigenvalues
        DeltaIv(:,:,i,f) = DeltaIv_; % insert the eigenvectors into array
        DeltaI(:,:,i,f) = DeltaI_; % insert the eigenvalues into array
    end
end

% -------------------------
% ---- gather results -----
% -------------------------
if options.statusBars
    waitbar(nframe/(nframe+1),h,'Gathering results....')
end
% general
m.x = x; % x1 x2 positions
m.FrameTimes = video.FrameTimes(video.Frames);
m.P = P; % raw data
% DDE
m.DDE.F = F; % DDE deformation gradient tensor
m.DDE.E = E; % DDE strain matrix
m.DDE.EIv = EIv; % DDE eigenvectors
m.DDE.EI = EI; % DDE eigenvalues
m.DDE.cond = Cond; % DDE matrix condition
% LSF
m.LSF.F = FLSF; % LSF deformation gradient tensor
m.LSF.E = ELSF;  % LSF strain matrix
m.LSF.EIv = ELSFIv;  % LSF strain eigenvectors matrix
m.LSF.EI = ELSFI;  % LSF strain eigenvalues matrix
m.LSF.cond = LSFCond; % LSF matrix condition
% SIMPLE
m.SIMPLE.Delta = Delta; % xy Delta
m.SIMPLE.DeltaIv = DeltaIv; % preallocate SIMPLE eigenvectors
m.SIMPLE.DeltaI = DeltaI; % preallocate SIMPLE eigenvalues
% averages
m.FrameAverage.F = FLSF_average;  % average deformation (LSF over all points)
m.FrameAverage.E = ELSF_average;  % average strain      (LSF over all points)
% save the results
mechanics = m;
% kill waitbar
if options.statusBars
    delete(h)
end
end

function  mechanics = processManifold2d(PC,video,P,X,options)
nd = size(X,2); % number of dimensions
% get the calibration parameters from the calibration matrix
K = calibrationParameterstoMatrix(options.CameraCalibration); 
% get the number of elements
nelement = size(P,1);
nframe = size(P,3);
% preallocate F, E, and X, FLSF, ELSF, and averages
F = zeros(nd,nd,nelement,nframe); % DDE deformation gradient tensor
E = zeros(nd,nd,nelement,nframe); % DDE strain matrix
x = zeros(nd,nelement,nframe); % x1 x2 positions
EIv = zeros(nd,nd,nelement,nframe); % DDE eigenvectors
EI = zeros(nd,nd,nelement,nframe); % DDE eigenvalues
%FLSF = zeros(nd,nd,nelement,nframe); % LSF deformation gradient tensor
%ELSF = zeros(nd,nd,nelement,nframe); % LSF strain matrix
%ELSFIv = zeros(nd,nd,nelement,nframe); % LSF strain eigenvectors matrix
%ELSFI = zeros(nd,nd,nelement,nframe); % LSF strain eigenvalues matrix
%FLSF_average = zeros(nd,nd,nframe); % average deformation (LSF over all points)
%ELSF_average = zeros(nd,nd,nframe); % average strain      (LSF over all points)
DeltaIv = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvectors
DeltaI = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvalues
% go through each element and each frame and calculate the deformation, the position, and
% the strain
if options.statusBars
    h = waitbar(0,['1 /' num2str(nframe)],'Name',['Processing Mechanics for ' num2str(nframe) ...
        ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
for f = 1:nframe
    if options.statusBars && getappdata(h,'canceling')% check if the cancel button has been pushed
        break % if so quit
    end
    % -------------------------
    % ----- DDE mechanics -----
    % -------------------------
    for i = 1:nelement
        Ft = feval(options.warpingFcn,P(i,:,f)); % use the Warping Function to generate a deformation gradient tensor for this Region
        xt = K.inv*PC(i).N.inv*Ft*PC(i).N.fwd*K.fwd*[X(i,:) 1]'; % warp the coordinates of X into x for LSF.
        Ft = Ft(1:nd,1:nd); % extract the deformation gradient Tensor (2d)
        F(:,:,i,f) = Ft; % add the deformation to the
        x(:,i,f) = xt(1:nd)'; % insert temporary position into position array
        E_ = 0.5*(Ft'*Ft-eye(nd)); % calculate DDE strain
        E(:,:,i,f) = E_; % insert strain into array
        [EIv_, EI_] = eig(E_); % get the eigenvectors and eigenvalues of the strain
        EIv(:,:,i,f) = EIv_; % insert the eigenvectors into array
        EI(:,:,i,f) = EI_; % insert the eigenvalues into array
    end
    % -------------------------
    % ----- LSF mechanics -----
    % -------------------------
%     xt = x(:,:,f)'; % Transpose x (to match X)
%     % find elements that are close to one and other to do LSF mechanics
%     LSF_groups = knnsearch(xt,xt,'K',options.mechanics.LSF_Points);
%     for i = 1:nelement
%         g = LSF_groups(i,:); % grab the 5 closest points
%         FLSF_ = LSF_F(X(g,:),xt(g,:)); % calculate the deformation
%         ELSF_ = 0.5*(FLSF_'*FLSF_-eye(nd)); % calculate strain based on deformation
%         FLSF(:,:,i,f) = FLSF_; % insert temp variable into entire array
%         ELSF(:,:,i,f) = ELSF_; % insert temp variable into entire array
%         [ELSFIv_,ELSFI_] = eig(ELSF_); % get LSF eigenvalues
%         ELSFIv(:,:,i,f) = ELSFIv_; % insert the eigenvectors into array
%         ELSFI(:,:,i,f) = ELSFI_; % insert the eigenvalues into array
%     end
%     % caclulate average strains based on displacements
%     FLSF_average_ = LSF_F(X,xt);
%     FLSF_average(:,:,f) = FLSF_average_;
%     ELSF_average(:,:,f) = 0.5*(FLSF_average_'*FLSF_average_-eye(nd));
    % update waitbar
    if options.statusBars
        waitbar(f/(nframe+1),h,[num2str(f) '/' num2str(nframe)])
    end
end

% -------------------------
% ---- gather results -----
% -------------------------
if options.statusBars 
    waitbar(nframe/(nframe+1),h,'Gathering results....')
end
% general
m.x = x; % x1 x2 positions
m.FrameTimes = video.FrameTimes(video.Frames);
% DDE
m.DDE.F = F; % DDE deformation gradient tensor
m.DDE.E = E; % DDE strain matrix
m.DDE.EIv = EIv; % DDE eigenvectors
m.DDE.EI = EI; % DDE eigenvalues
% LSF
%m.LSF.F = FLSF; % LSF deformation gradient tensor
%m.LSF.E = ELSF;  % LSF strain matrix
%m.LSF.EIv = ELSFIv;  % LSF strain eigenvectors matrix
%m.LSF.EI = ELSFI;  % LSF strain eigenvalues matrix
% SIMPLE
%m.SIMPLE.Delta = Delta; % xy Delta
%m.SIMPLE.DeltaIv = DeltaIv; % preallocate SIMPLE eigenvectors
%m.SIMPLE.DeltaI = DeltaI; % preallocate SIMPLE eigenvalues
% averages
%m.FrameAverage.F = FLSF_average;  % average deformation (LSF over all points)
%m.FrameAverage.E = ELSF_average;  % average strain      (LSF over all points)
% save the results
mechanics = m;
end

function mechanics = processAffine1d(PC,video,P,X,options)
nd = size(X,2); % number of dimensions
% get the calibration parameters from the calibration matrix
K = calibrationParameterstoMatrix(options.CameraCalibration); 
% get the number of elements
nelement = size(P,1);
nframe = size(P,3);
% preallocate F, E, and X, FLSF, ELSF, and averages
F = zeros(nd,nd,nelement,nframe); % DDE deformation gradient tensor
E = zeros(nd,nd,nelement,nframe); % DDE strain matrix
x = zeros(nd,nelement,nframe); % x1 x2 positions
EIv = zeros(nd,nd,nelement,nframe); % DDE eigenvectors
EI = zeros(nd,nd,nelement,nframe); % DDE eigenvalues
FLSF = zeros(nd,nd,nelement,nframe); % LSF deformation gradient tensor
ELSF = zeros(nd,nd,nelement,nframe); % LSF strain matrix
ELSFIv = zeros(nd,nd,nelement,nframe); % LSF strain eigenvectors matrix
ELSFI = zeros(nd,nd,nelement,nframe); % LSF strain eigenvalues matrix
FLSF_average = zeros(nd,nd,nframe); % average deformation (LSF over all points)
ELSF_average = zeros(nd,nd,nframe); % average strain      (LSF over all points)
DeltaIv = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvectors
DeltaI = zeros(nd,nd,nelement,nframe); % preallocate SIMPLE eigenvalues
Cond = zeros(nelement,nframe);
LSFCond = zeros(nelement,nframe);
% go through each element and each frame and calculate the deformation, the position, and
% the strain
if options.statusBars
    h = waitbar(0,['1 /' num2str(nframe)],'Name',['Processing Mechanics for ' num2str(nframe) ...
        ' frames...'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
end
for f = 1:nframe
    if options.statusBars && getappdata(h,'canceling')% check if the cancel button has been pushed
        break % if so quit
    end
    % -------------------------
    % ----- DDE mechanics -----
    % -------------------------
    for i = 1:nelement
        Ft = feval(options.warpingFcn,P(i,:,f)); % use the Warping Function to generate a deformation gradient tensor for this Region
        xt = Ft*PC(i).N.fwd*K.fwd*[X(i,:) 1]'; % warp the coordinates of X into x for LSF.
        xt = xt./xt(end);
        xt = K.inv*PC(i).N.inv*xt;
        Ft = Ft(1:nd,1:nd); % extract the deformation gradient Tensor (2d)
        F(:,:,i,f) = Ft; % add the deformation to the
        x(:,i,f) = xt(1:nd)./repmat(xt(nd+1),nd,1); % insert temporary position into position array
        E_ = 0.5*(Ft'*Ft-eye(nd)); % calculate DDE strain
        E(:,:,i,f) = E_; % insert strain into array
        [EIv_, EI_] = eig(E_); % get the eigenvectors and eigenvalues of the strain
        EIv(:,:,i,f) = EIv_; % insert the eigenvectors into array
        EI(:,:,i,f) = EI_; % insert the eigenvalues into array
        Cond(i,f) = rcond(Ft); % get condition matrix for DDE
    end
    % -------------------------
    % ----- LSF mechanics -----
    % -------------------------
    xt = x(:,:,f)'; % Transpose x (to match X)
    % find elements that are close to one and other to do LSF mechanics
    LSF_groups = [1:nelement-1; 2:nelement]';
    for i = 1:nelement-1
        g = LSF_groups(i,:); % grab the closest points
        Xd = hypot(X(g(1),1)-X(g(2),1),X(g(1),2)-X(g(2),2));
        xd = hypot(xt(g(1),1)-xt(g(2),1),xt(g(1),2)-xt(g(2),2));
        xlam(:,i,f) = [xt(g(1),1)+xt(g(2),1) xt(g(1),2)+xt(g(2),2)]/2;
        lam(i,f) = xd/Xd;
        Elam(i,f)  = 0.5*(lam(i,f)^2-1); % calculate strain based on deformation
    end
    % caclulate average strains based on displacements
    FLSF_average_ = LSF_F(X,xt);
    FLSF_average(:,:,f) = FLSF_average_;
    ELSF_average(:,:,f) = 0.5*(FLSF_average_'*FLSF_average_-eye(nd));
    % update waitbar
    if options.statusBars
        waitbar(f/(nframe+1),h,[num2str(f) '/' num2str(nframe)])
    end
end


% -------------------------
% ---- gather results -----
% -------------------------
if options.statusBars
    waitbar(nframe/(nframe+1),h,'Gathering results....')
end
% general
m.x = x; % x1 x2 positions
m.FrameTimes = video.FrameTimes(video.Frames);
m.P = P; % raw data
% DDE
m.DDE.F = F; % DDE deformation gradient tensor
m.DDE.E = E; % DDE strain matrix
m.DDE.EIv = EIv; % DDE eigenvectors
m.DDE.EI = EI; % DDE eigenvalues
m.DDE.cond = Cond; % DDE matrix condition

% LSF
m.LSF.x = xlam;
m.LSF.F = lam; % LSF deformation gradient tensor
m.LSF.E = Elam;  % LSF strain matrix

% averages
m.FrameAverage.F = FLSF_average;  % average deformation (LSF over all points)
m.FrameAverage.E = ELSF_average;  % average strain      (LSF over all points)
% save the results
mechanics = m;
% kill waitbar
if options.statusBars
    delete(h)
end
end

function F = LSF_F(X,x)
% calculate deformation gradient tensor based on the formula
% dr = F*dR
% X = points at undeformed configuration
% x = points at deformed configuration

dR = bsxfun(@minus,X,mean(X,1));
dr = bsxfun(@minus,x,mean(x,1));
F = dR\dr;
end



