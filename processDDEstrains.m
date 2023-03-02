%% Change strains 

%% Saving .csv file *may need to re-grid due to transformation done by DDE
testName = input('What do you want to name this test?: ','s');
saveDir = input('What directory would you like to save to?: ', 's');
DDEcode = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\DDECode';

% Reference Configuration (same for every frame, outside loop)
x = mechanics.x( 2, :, 1 )';
y = mechanics.x( 1, :, 1 )';

dims = size( mechanics.x );
dataPoints = dims( 2 ); % number of pixels evaluated
nframes = dims( 3 ); % number of frames evaluated
sigma = zeros( 1, dataPoints )';
gamma = zeros( 1, dataPoints )';


cd( saveDir )
save( 'mechanics.mat', 'mechanics' );
wb = waitbar( 0, 'Starting DDE .csv writing...' );
for i = 1:nframes
    waitbar( i/nframes, wb, ...
        'Please wait for DDE .csv writing...' );
    fileName = sprintf( '%s%04d', [testName '_'], i ); %name for .csv file
    
    % Current Configuration
    x_c = mechanics.x( 2, :, i )';
    y_c = mechanics.x( 1, :, i )';
    
    % Displacements
    u = x_c - x; % x-displacement
    v = y_c - y; % y-displacement
    u_c = u; v_c = v; % not sure what these are, won't use
    
    % Strains
    exx = zeros( 1, dataPoints )'; % longitudinal strain vector
    eyy = zeros( 1, dataPoints )'; % transverse strain vector
    exy = zeros( 1, dataPoints )'; % shear strain vector
    e1 = zeros( 1, dataPoints )'; % maximum principal strain vector
    e2 = zeros( 1, dataPoints )'; % minimum principal strain vector
    for j = 1:dataPoints
        defGrad = mechanics.DDE.F( :, :, j, i ); % Deformation Gradient Tensor
        strainGL = 0.5*(defGrad'*defGrad - eye( 2, 2 )); %Green-Lagrange Strain Tensor (a.k.a Material Strain Tensor)
        exx( j, 1 ) = strainGL( 2, 2 );
        eyy( j, 1 ) = strainGL( 1, 1 );
        exy( j, 1 ) = strainGL( 1, 2 );
    end % principal strains are calculated in DispErrorCalc code
    
    % need to create the file
    T = table( x, y, u, v, sigma, x_c, y_c, u_c, v_c, exx, eyy, exy, e1, e2, gamma,...
        'VariableNames',{'x', 'y', 'u', 'v', 'sigma', 'x_c', 'y_c', 'u_c', 'v_c', 'exx', 'eyy', 'exy', 'e1', 'e2', 'gamma'} );
    writetable(T, [fileName '.csv']);
end
close( wb )
cd( DDEcode )