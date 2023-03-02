%% Main DDE Code, adapted from Boyle et al.
%% Adaptations by Frederick Houghton (houghton.frederick@gmail.com)
%% Last Updated 06/07/22

inputType = '3dmat'; % setup the input type
% folder containing all of the mat files of hearts
% directory = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\GUI Results\Noiseless'; 
% directory = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study';
directory = '/Users/fhoughton/Desktop/DDE Final Code/Example';
getFrameTimes = false; % no frame times for these videos
warning off % disable warnings
% files to process
files = {'10DEF_Noiseless_dataIMG'};

% Hard Coded Variables
targetSize = [400 200 300]; % size of mat file
subsetDia = 81; %assuming isotropic search
stepSize = 14; %assuming isotropic search

% designing AOI based on subset size: to borders of image
imHeight = targetSize( 1 );
imWidth = targetSize( 2 );
boxRad = (subsetDia-1)/2; % calculate distance from node to edge of subset
rImHeight = imHeight - 2*boxRad;
rImWidth = imWidth - 2*boxRad;

hAOI = rImHeight/imHeight; wAOI = rImWidth/imWidth;


options = setLKoptions('WarpType'                           ,'affine2d', ... 'affine3d'
                       'ImageProcessing.InterpolantWarnings',false, ... 
                       'Meshing.boxSize'                    ,[subsetDia subsetDia], ...  [11 11]                     
                       'Meshing.spacing'                    ,[stepSize stepSize], ... [6 6]
                       'Optimization.maximumIterations'     , 80 , ...
                       'CameraCalibration'                  ,'off', ...
                       'ImageProcessing.demosaic.on'        ,false, ...
                       'Optimization.epsilon'               ,.006 , ...
                       'Meshing.AnalysisRegionPrecentage'   ,[hAOI wAOI] , ... [0.75 0.75]
                       'statusBars'                         , false, ...save
                       'Meshing.mode'                       ,'centerBoxes', ...
                       'mechanics.LSF_Points'               ,7        , ...  
                       'ImageProcessing.FrameResize.Ratio'  , 1     , ...
                       'mechanics.fixDim'                   , false)         ; 
%% analyze each video
for ind = 1:numel(files) %loop through each heart
    tic % time it 
    file = [directory '/' files{ind} '.mat']; % load the mat file
    video = initializeVideo(inputType,getFrameTimes,file ); % intialize the video
    resizeRat = (prod(targetSize)/prod(video.size(1:3)))^(1/3); % adjust the numer of voxels to be approximately the same size as target
    options.ImageProcessing.FrameResize.Ratio = resizeRat; % resize the video to get same number of voxels as target 
    mechanics = lkInverseCompositional(video,options); % analyze the heart 
    tim = toc; % record the time
    fprintf('Sample %s processed in %dm\n',files{ind},round(tim/60)) % display output
    tic % time how long it takes to save the results 
    save(['allMatch' files{ind}]) % save the results
    tim2 = toc; % record the time
    fprintf('Sample %s saved in %dm\n',files{ind},round(tim2/60)) % print the time
end
%% mesh heatmaps and plot infarctions
% % specify files to analyze
% files = {'10DEF_Noiseless_dataIMG'};
% 
% % for each sample mesh the data and plot the result
% for ind = 1:numel(files)
%     % 0. load the data
%     fprintf('loading data for %s\n',files{ind}) % display info 
%     load(['allMatch' files{ind}],'mechanics','video'); % load the actual data
%     % 1. mesh th resulting data
%     clear Vh c % clear these to avoid problems
%     fprintf('generating mesh for %s\n',files{ind})
%     Vh = generateHeartHeatmap_2d_interp3d(video,mechanics); % this step is long
%     keyboard;
%     % 2. generate an output plot
%     fprintf('plotting results for %s\n',files{ind})
%     c = infarcPlottingLoop(files{ind},video,mechanics.options,Vh,curDat);
% end

%% Plot strains
% close all 
% frame_no = 500;
% x = mechanics.x(1,:,1)';
% y = mechanics.x(2,:,1)';
% exx = squeeze(mechanics.DDE.E(2,2,:,frame_no));
% eyy = mechanics.DDE.E(1,1,:,frame_no);
% exy = mechanics.DDE.E(1,2,:,frame_no);
% 
% strain = load('C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\Strains_XX\frame0500_strain_xx.mat');
% strain = strain.strain;
% 
% figure
% subplot(1,2,1)
% scatter(y,x,500,exx,'filled','s');
% set(gca,'YDir','reverse') 
% h = video.size(1);
% w = video.size(2);
% pbaspect([w,h,1]);
% colorbar
% caxis([0 0.1])
% xlim([1 w])
% ylim([1 h])
% 
% 
% subplot(1,2,2)
% imagesc(strain)
% pbaspect([w,h,1]);
% colorbar
% caxis([0 0.1])