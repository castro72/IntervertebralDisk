%% Auto-Run DDE Code
%% By Frederick Houghton (houghton.frederick@gmail.com)
%% Last Updated 06/07/22

% toolbox = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Toolbox';
toolbox = '/Users/fhoughton/Desktop/Benchmark Toolbox';
% DDEcode = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\DDECode';
DDEcode = '/Users/fhoughton/Desktop/DDEcode';

%% Instantiate
targetSize = [400 200 1101];
subsets = [81];
    % [11 21 31 41 51 61 71 81 91 101]
stepSizes = [14];
    % [2 6 10 14];
    
%% Noiseless 10DEF    
directory = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF';
files = {'10DEF_Noiseless_dataIMG'}; %#ok<*NASGU>
for ij = 1:length( subsets )
    for jk = 1:length( stepSizes )
        subsetDia = subsets( ij ) %#ok<*NOPTS>
        stepSize = stepSizes( jk )
        
        cd( DDEcode )
        run('analyzeAllHeartsMatchedSize_2D_FreddieAndReece.m')
        
        compTimefolder = sprintf( '%s%d%s%d', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noiseless\SS', subsetDia, '_S', stepSize );
        cd( compTimefolder )
        tim = tim/60; % convert to minutes
        save( 'Computational Effort.mat', 'tim' )
        testName = sprintf( '%s%d%s%d', '10DEF_NN_SS', subsetDia, '_S', stepSize );
        saveDir = sprintf( '%s%d%s%d%s', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noiseless\SS', subsetDia, '_S', stepSize, '\CSV Files' );
        cd( DDEcode )
        run('processDDEstrains.m')
    end
end

%% Noise 10DEF
% directory = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF';
directory = '/Volumes/Tendon Strain/Tendon Strain Study/Benchmark Results/Final Synthetic Test Cases/10DEF';
files = {'10DEF_Noise_dataIMG'};
for ij = 1:length( subsets )
    for jk = 1:length( stepSizes )
        subsetDia = subsets( ij )
        stepSize = stepSizes( jk )
        
        cd( DDEcode )
        run('analyzeAllHeartsMatchedSize_2D_FreddieAndReece.m')
        
%         compTimefolder = sprintf( '%s%d%s%d', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noise\SS', subsetDia, '_S', stepSize );
        compTimefolder = sprintf( '%s%d%s%d', '/Volumes/Tendon Strain/Tendon Strain Study/Benchmark Results/Final Synthetic Test Cases/10DEF/DDE Results/Noise/SS', subsetDia, '_S', stepSize );
        tim = tim/60; % convert to minutes
        save( 'Computational Effort.mat', 'tim' )
        testName = sprintf( '%s%d%s%d', '10DEF_YN_SS', subsetDia, '_S', stepSize );
%         saveDir = sprintf( '%s%d%s%d%s', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noise\SS', subsetDia, '_S', stepSize, '\CSV Files' );
        saveDir = sprintf( '%s%d%s%d%s', '/Volumes/Tendon Strain/Tendon Strain Study/Benchmark Results/Final Synthetic Test Cases/10DEF/DDE Results/Noise/SS', subsetDia, '_S', stepSize, '/CSV Files' );
        run('processDDEstrains.m')
    end
end

%% Noiseless 10DEF Error Calc DDE
cd( toolbox )
for ij = 1:length( subsets )
    for jk = 1:length( stepSizes )
        subsetDia = subsets( ij )
        stepSize = stepSizes( jk )
        DDE_folder = sprintf( '%s%d%s%d%s', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noiseless\SS', subsetDia, '_S', stepSize, '\CSV Files' );
        save_folder = sprintf( '%s%d%s%d', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noiseless\SS', subsetDia, '_S', stepSize' );
        strain_folder = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF';
        frame_rate = 100;
        
        run('DDE_DispErrorCalc.m')
        clearvars -except toolbox DDEcode directory files targetSize subsets...
            stepSizes ij jk
    end
end

%% Noise 10DEF Error Calc DDE
cd( toolbox )
for ij = 1:length( subsets )
    for jk = 1:length( stepSizes )
        subsetDia = subsets( ij )
        stepSize = stepSizes( jk )
        DDE_folder = sprintf( '%s%d%s%d%s', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noise\SS', subsetDia, '_S', stepSize, '\CSV Files' );
        save_folder = sprintf( '%s%d%s%d', 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF\DDE Results\Noise\SS', subsetDia, '_S', stepSize' );
        strain_folder = 'C:\Users\gocon\OneDrive\Desktop\Tendon Strain Study\Benchmark Results\Final Synthetic Test Cases\10DEF';
        frame_rate = 100;
        
        run('DDE_DispErrorCalc.m')
        clearvars -except toolbox DDEcode directory files targetSize...
            subsets stepSizes ij jk
    end
end