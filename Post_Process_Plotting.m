%% Post-Processing Plots
%% by Frederick Houghton (houghton.frederick@gmail.com)
%% Last Updated: June 9, 2022

% direct parsing of data received from DDE code
DDE_folder = input( 'What folder are your .csv files in?: ', 's' );
save_folder = input( 'Where would you like to save the results of this code?: ', 's' );
step_size = input( 'What step size was used in the DDE code for this test?: ' );
[x, y, exx, eyy, exy, e1, e2] = DDEParse( DDE_folder, step_size );

% read-in the .mat file that DDE was performed on
dataFile = input( 'Please input the location and name of the data file you would like to use:', 's' );
load( dataFile ); %you will have to hard code the name of your data file into this code

%% Show Video or Build .avi files (can change variables if desired)
choice = input('Would you like to build the .avi file? (Y/N) : ', 's' );
if choice == 'Y'
    %threshold for color scales in images
    thresh = [0 0.1].*110;
    frame_count = 300;
    % Exx video
    exxMap = VideoWriter( sprintf( '%s%s', ['/Users/fhoughton/Desktop','/'],...
        'exxMap.avi' ) );
%     exxMap = VideoWriter( sprintf( '%s%s', [save_folder,'/'],...
%         'exxMap.avi' ) );
    exxMap.FrameRate = frame_count/60;
    open( exxMap );
    path = save_folder;
    for i = 1 : frame_count
        f1 = figure(1);
        %raw data print out
        p1 = pcolor( exx( :, :, i ).*100 );
        p1.EdgeColor = 'none';
        colormap( 'default' )
        caxis( thresh )
        colorbar
        sgtitle( sprintf( '%s%04d', 'Frame Number: ', i ) )
        frame = getframe( gcf );
        writeVideo( exxMap, frame );
    end
    close( exxMap ); clear frame
    
    % Eyy video
%     eyyMap = VideoWriter( sprintf( '%s%s', [save_folder,'\'],...
%         'eyyMap.avi' ) );
    eyyMap = VideoWriter( sprintf( '%s%s', ['/Users/fhoughton/Desktop','/'],...
        'eyyMap.avi' ) );
    eyyMap.FrameRate = frame_count/60;
    open( eyyMap );
    for j = 1 : frame_count
        f2 = figure(2);
        p1 = pcolor( DIC_eyy( :, :, j ).*100 );
        p1.EdgeColor = 'none';
        colormap( 'default' )
        caxis( thresh )
        colorbar
        sgtitle( sprintf( '%s%04d', 'Frame Number: ', j ) )
        frame = getframe( gcf );
        writeVideo( eyyMap, frame );
    end
    close( eyyMap );
end
clear i j

%% Strain Maps and Data Overlay
frame = input( 'What frame would you like to plot?: ' );

exx_frame = exx( :, :, frame );
eyy_frame = eyy( :, :, frame );
exy_frame = exy( :, :, frame );

x_range = [min( x ) max( x )]; 
y_range = [min( y ) max( y )]; 
x = x_range( 1 ):step_size:x_range( 2 );
y = y_range( 1 ):step_size:y_range( 2 );

transparency = zeros( size( dataIMG( :, :, frame ) ) );
cols = x_range( 1 ):x_range( 2 ); rows = y_range( 1 ):y_range( 2 );
transparency( rows, cols ) = 1;

% interpolating strains for plots
[Xq,Yq] = meshgrid( cols, rows );
exx_interp = interp2( x, y, exx_frame, Xq, Yq );
eyy_interp = interp2( x, y, eyy_frame, Xq, Yq );
exy_interp = interp2( x, y, exy_frame, Xq, Yq );

% add NaN values to borders to fit image over data
exx_plot = nan( size( dataIMG( :, :, frame ) ) );
eyy_plot = nan( size( dataIMG( :, :, frame ) ) );
exy_plot = nan( size( dataIMG( :, :, frame ) ) );
exx_plot( rows, cols ) = exx_interp;
eyy_plot( rows, cols ) = eyy_interp;
exy_plot( rows, cols ) = exy_interp;
opacity = 0.5; %hard code to your liking

BWframe = repmat( uint8( dataIMG( :, :, frame ) ), [1 1 3] ); %expand BW image

f3 = figure( 3 ); %Exx strain map
imshow( BWframe ); hold on
hexx_plot = imshow( exx_plot );
hexx_plot.AlphaData = transparency*opacity;
colormap( hsv )
caxis( [-0.15 0.15] );
colorbar;

f4 = figure( 4 ); %Eyy strain map
imshow( BWframe ); hold on
heyy_plot = imshow( eyy_plot );
heyy_plot.AlphaData = transparency*opacity;
colormap( hsv )
caxis( [-0.15 0.15] );
colorbar;

f5 = figure( 5 ); %Exy strain map
imshow( BWframe ); hold on
hexy_plot = imshow( exy_plot );
hexy_plot.AlphaData = transparency*opacity;
colormap( hsv )
caxis( [-0.15 0.15] );
colorbar;

%% Standard Plots
[rows, cols, num_frames] = size( exx );
% Medians and means and poisson's ratios
med_exx = zeros( num_frames, 1 ); med_eyy = zeros( num_frames, 1 );
med_exy = zeros( num_frames, 1 ); mean_exx = zeros( num_frames, 1 );
mean_eyy = zeros( num_frames, 1 ); mean_exy = zeros( num_frames, 1 );
poisson_mean = zeros( num_frames, 1 );poisson_med = zeros( num_frames, 1 );
for i = 1:num_frames
    med_exx( i, 1 ) = median( median( exx( :, :, i ), 'omitnan' ) );
    med_eyy( i, 1 ) = median( median( eyy( :, :, i ), 'omitnan' ) );
    med_exy( i, 1 ) = median( median( exy( :, :, i ), 'omitnan' ) );
    mean_exx( i, 1 ) = mean( mean( exx( :, :, i ), 'omitnan' ) );
    mean_eyy( i, 1 ) = mean( mean( eyy( :, :, i ), 'omitnan' ) );
    mean_exy( i, 1 ) = mean( mean( exy( :, :, i ), 'omitnan' ) );
    poisson_mean( i, 1 ) = -mean_eyy( i, 1 )/mean_exx( i, 1 );
    poisson_med( i, 1 ) = -med_eyy( i, 1 )/med_exx( i, 1 );
end

frameVec = linspace( 1, 300, 300 )';
med_e1 = median( e1, 1, 'omitnan' )'; mean_e1 = mean( e1, 1, 'omitnan' )';
med_e2 = median( e2, 1, 'omitnan' )'; mean_e2 = mean( e2, 1, 'omitnan' )';

prin_poisson_med = -med_e2./med_e1;
prin_poisson_mean = -mean_e2./mean_e1;

f8 = figure( 8 ); % Exx and Eyy over time
plot( frameVec, med_exx, 'ro' )
hold on
plot( frameVec, mean_exx, 'bo' )
plot( frameVec, med_eyy, 'o' )
plot( frameVec, mean_eyy, 'o' )
xlabel( 'Frame Number', 'Interpreter', 'latex', 'FontSize', 16 )
ylabel( 'Strain', 'Interpreter', 'latex', 'FontSize', 16 )
title( 'Longitudinal and Transverse Strains', 'Interpreter', 'latex', 'FontSize', 24 )
legend( 'Median Exx', 'Mean Exx', 'Median Eyy', 'Mean Eyy', 'Interpreter', 'latex', 'FontSize', 16 )

f9 = figure( 9 ); % Exy over time
plot( frameVec, med_exy, 'ro' )
hold on
plot( frameVec, mean_exy, 'bo' )
xlabel( 'Frame Number', 'Interpreter', 'latex', 'FontSize', 16 )
ylabel( 'Strain', 'Interpreter', 'latex', 'FontSize', 16 )
title( 'Shear Strain', 'Interpreter', 'latex', 'FontSize', 24 )
legend( 'Median', 'Mean', 'Interpreter', 'latex', 'FontSize', 16 )

f10 = figure( 10 ); % Principal Strains over time
plot( frameVec, med_e1, 'ro' )
hold on
plot( frameVec, mean_e1, 'bo' )
plot( frameVec, med_e2, 'go' )
plot( frameVec, mean_e2, 'o' )
xlabel( 'Frame Number', 'Interpreter', 'latex', 'FontSize', 16 )
ylabel( 'Strain', 'Interpreter', 'latex', 'FontSize', 16 )
title( 'Principal Strains', 'Interpreter', 'latex', 'FontSize', 24 )
legend( 'Median E1', 'Mean E1', 'Median E2', 'Mean E2', 'Interpreter', 'latex', 'FontSize', 16 )

f11 = figure( 11 ); %Poisson's Ratio plots
plot( frameVec, poisson_med, 'o' )
hold on
plot( frameVec, poisson_mean, 'o' )
plot( frameVec, prin_poisson_med, 'o' )
plot( frameVec, prin_poisson_mean, 'o' )
xlabel( 'Frame Number', 'Interpreter', 'latex', 'FontSize', 16 )
ylabel( 'Ratio', 'Interpreter', 'latex', 'FontSize', 16 )
title( 'Poissons Ratio', 'Interpreter', 'latex', 'FontSize', 24 )
legend( 'Median $\frac{-Eyy}{Exx}$', 'Mean $\frac{-Eyy}{Exx}$', ...
    'Median $\frac{-E2}{E1}$', 'Mean $\frac{-E2}{E1}$', 'Interpreter',...
    'latex', 'FontSize', 16 )

%% Parsing Function
function [x, y, DDE_exx, DDE_eyy, DDE_exy, DDE_e1, DDE_e2] = DDEParse(DDE_folder,step_size) 
% DDE_folder = input( 'Choose the directory to use for DDE2D: ', 's' );
% save_folder = input( 'Choose a directory to save to: ', 's' );
% strain_folder = input( 'Choose a directory to locate true values: ', 's' );
% frame_rate = input( 'Frame rate?: ' );
% step_size = input( 'Step Size?: ' );

start_folder = cd;
cd( DDE_folder )
fprintf( 'Reading files from %s\n', pwd )
fileList = dir('*.csv');

% Finding the size of the mapping using reference case
reffile = fileList( 1 ).name; refdata = csvread( reffile, 1, 0 );
% Note: DDE2D orders data from left to right then top-down
x = refdata( :, 1 ); y = refdata( :, 2 ); % get locations
coords = [x y]; % store for future use
% get size of matrix
cols = ((x( end ) - x( 1 ))/step_size) + 1;
rows = ((y( end ) - y( 1 ))/step_size) + 1;
frame_count = length( fileList );

u = nan( rows, cols, length( fileList ) );
v = nan( rows, cols, length( fileList ) );
exx = nan( rows, cols, length( fileList ) );
eyy = nan( rows, cols, length( fileList ) );
exy = nan( rows, cols, length( fileList ) );
DDE_u = nan( rows, cols, length( fileList ) );
DDE_v = nan( rows, cols, length( fileList ) );
DDE_exx = nan( rows, cols, length( fileList ) );
DDE_eyy = nan( rows, cols, length( fileList ) );
DDE_exy = nan( rows, cols, length( fileList ) );
exx_flag = zeros( rows*cols, length(fileList ) );
eyy_flag = zeros( rows*cols, length(fileList ) );
exy_flag = zeros( rows*cols, length(fileList ) );
clear refdata reffile

wb = waitbar( 0, 'Starting DDE parsing...' );

for i = 1:length( fileList )
    waitbar( i/length( fileList ), wb, ...
        'Please wait for DDE output parsing...' );
    fileName = fileList( i ).name;
    data = []; %#ok<*NASGU>
    try
        data = csvread( fileName, 1, 0 );
    catch
        fprintf( 'File %s is empty\n', fileName )
        
        % No value added to matrices at this index
        u( :, :, i ) = nan;
        v( :, :, i ) = nan;
        exx( :, :, i ) = nan;
        eyy( :, :, i ) = nan;
        exy( :, :, i ) = nan;
        %         e1( i ) = nan;
        %         e2( i ) = nan;
        %         gamma( i ) = nan;
        
        continue
    end
    
    % pull out displacement and strain vectors and store them as a map
    % for each frame
    u_frame = data( :, 3 ); v_frame = data( :, 4 );
    exx_frame = data( :, 10 ); eyy_frame = data( :, 11 );
    exy_frame = data( :, 12 );
    % e1 = data( :, 13 );
    % e2 = data( :, 14 );
    % gamma = data( :, 15 );
    
    % take absolute value of strains for filtering
    %     abs_u = abs( u_frame( : ) );
    %     abs_v = abs( v_frame( : ) );
    abs_exx = abs( exx_frame( : ) );
    abs_eyy = abs( eyy_frame( : ) );
    abs_exy = abs( exy_frame( : ) );
    
    s_max = 0.5; % threshold to filter strain values
    
    % find all absolute values exceeding 500% strain: NaN (interpolate?)
    %     u_ind = find( abs_u >= 5 ); u_frame( u_ind ) = NaN;
    %     v_ind = find( abs_v >= 5 ); v_frame( v_ind ) = NaN;
    exx_ind = find( abs_exx >= s_max ); exx_frame( exx_ind ) = NaN; %#ok<*FNDSB>
    eyy_ind = find( abs_eyy >= s_max ); eyy_frame( eyy_ind ) = NaN;
    exy_ind = find( abs_exy >= s_max ); exy_frame( exy_ind ) = NaN;
    
    % flag any values between 5 and 100
    %     u_flag = find( abs_u >= 5 && abs_u <= 100 );
    %     v_flag = find( abs_v >= 5 && abs_v <= 100 );
    for n = 1 : numel( exx_frame )
        if abs_exx( n ) >= s_max && abs_exx( n ) <= 100
            exx_flag( n, i ) = 1;
        end
        if abs_eyy( n ) >= s_max && abs_eyy( n ) <= 100
            eyy_flag( n, i ) = 1;
        end
        if abs_exy( n ) >= s_max && abs_exy( n ) <= 100
            exy_flag( n, i ) = 1;
        end
    end
    %     for n = 1 : rows*cols
    %         if abs_exx( n ) >= s_max
    %             exx_flag( n, i ) = 1;
    %         end
    %         if abs_eyy( n ) >= s_max
    %             eyy_flag( n, i ) = 1;
    %         end
    %         if abs_exy( n ) >= s_max
    %             exy_flag( n, i ) = 1;
    %         end
    %     end
    
    % write vectors for i-th frame to the arrays for each variable
    count = 1;
    for j = 1 : rows
        for k = 1 : cols
            u( j, k, i ) = u_frame( count, 1 );
            v( j, k, i ) = v_frame( count, 1 );
            exx( j, k, i ) = exx_frame( count, 1 );
            eyy( j, k, i ) = eyy_frame( count, 1 );
            exy( j, k, i ) = exy_frame( count, 1 );
            count = count + 1;
        end
    end
end
close( wb )
cd( start_folder )


% filter out any Inf values (NaN is default for each array: not filtered)
for k = 1 : frame_count
    %     exx_cutoff = exx(:,:,k);
    %     exx_cutoff = exx_cutoff(:);
    for j = 1 : cols
        for i = 1 : rows
            if isfinite( u( i, j, k ) )
                DDE_u( i, j, k ) = u( i, j, k );
            end
            if isfinite( v( i, j, k ) )
                DDE_v( i, j, k ) = v( i, j, k );
            end
            if isfinite( exx( i, j, k ) )
                DDE_exx( i, j, k ) = exx( i, j, k );
            end
            if isfinite( eyy( i, j, k ) )
                DDE_eyy( i, j, k ) = eyy( i, j, k );
            end
            if isfinite( exy( i, j, k ) )
                DDE_exy( i, j, k ) = exy( i, j, k );
            end
        end
    end
end

% store displacements as one array
DDE_displacements = [];
DDE_displacements( :, :, :, 1 ) = DDE_u;
DDE_displacements( :, :, :, 2 ) = DDE_v;


DDE_e1 = zeros( rows*cols, frame_count );
DDE_e2 = zeros( rows*cols, frame_count );

wb = waitbar( 0, 'Please wait for principal strain calculations...' );
for kk = 1 : frame_count
    count = 1;
    for ii = 1 : rows
        for jj = 1 : cols
            % strain tensor for given point
            DDE_strainTensor = [DDE_exx( ii, jj, kk ), ...
                DDE_exy( ii, jj, kk ); DDE_exy( ii, jj, kk ),...
                DDE_eyy( ii, jj, kk )];
            
            % Principal Strains
            DDE_e1( count, kk ) = (DDE_strainTensor( 1, 1 ) +...
                DDE_strainTensor( 2, 2 ))/2 + sqrt( ...
                ((DDE_strainTensor( 1, 1 ) - DDE_strainTensor( 2, ...
                2 ))/2)^2 + DDE_strainTensor( 1, 2 )^2 );
            DDE_e2( count, kk ) = (DDE_strainTensor( 1, 1 ) +...
                DDE_strainTensor( 2, 2 ))/2 - sqrt( ...
                ((DDE_strainTensor( 1, 1 ) - DDE_strainTensor( 2, ...
                2 ))/2)^2 + DDE_strainTensor( 1, 2 )^2 );
            
            count = count + 1;
        end
    end
    waitbar( kk/length( fileList ), wb, ...
        'Please wait for principal strain calculations...' );
end
close( wb )

end